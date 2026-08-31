#!/usr/bin/env python
# ----------------------------------------------------------
# Time-series quality control: realistic-sigma error scaling and
# iterative outlier editing.
#
# Ports two mechanisms from GAMIT/GLOBK ``tsfit``:
#   * ``real_stats``  (gen_util/real_stats.f)  -> realistic_sigma
#   * ``edit_ns``     (tsfit.f)                 -> flag_outliers
#
# All functions are pure and stateless so they can be unit-tested in
# isolation. The tsfitting engine wires them together in doFitting().
# ----------------------------------------------------------
import logging
import numpy as np

# Default fitting/quality-control options. Kept as the single source of
# defaults when a caller does not pass ``fit_opts`` explicitly.
DEFAULT_FIT_OPTS = {
    # Uncertainty scaling of the covariance returned by curve_fit:
    #   'none'      -- use the raw N^-1 (absolute_sigma=True), no scaling;
    #   'nrms'      -- the current behaviour (scipy default, == N^-1 * chi2/dof);
    #   'realistic' -- Herring bin-averaging + FOGM extrapolation (real_stats).
    'sigma_scale'   : 'nrms',

    # realistic-sigma guards (mirrors tsfit's MINNUM).
    'min_rsig'      : 30,    # minimum number of used points to attempt FOGM
    'min_sigscale'  : 10,    # minimum dof below which sig_scale is forced to 1

    # Iterative outlier editing (mirrors tsfit's NSIGMA / edit_ns).
    'outlier'       : False,  # enable post-fit n-sigma editing
    'nsigma'        : 4.0,    # rejection threshold
    'outlier_scale' : 'mad',  # 'mad' (robust) or 'sigma' (tsfit's sigma*nrms)
    'max_iter'      : 10,     # maximum edit/refit passes
    'restore_factor': 0.9,    # hysteresis band for restoring points

    # Pre-fit sigma screening (mirrors tsfit's MAX_SIGMA, in mm here).
    # None / 0.0 disables the screen.
    'max_sigma'     : None,
}

# Trial FOGM correlation times scanned by real_stats (days).
_TAUS = np.array([1, 2, 4, 8, 16, 32, 64, 128, 256, 512,
                  1024, 2048, 4096, 8192, 16384], dtype=float)

# Upper bound on the number of averaging intervals (real_stats's max_av).
_MAX_AV = 150


def merge_fit_opts(opts=None):
    '''
    Merge caller-supplied options over the defaults, returning a fresh dict.

    Input:
        opts = dict of fit/quality-control options, or None for all defaults

    Output:
        a dict containing every key of DEFAULT_FIT_OPTS
    '''
    merged = dict(DEFAULT_FIT_OPTS)
    if opts:
        for key, value in opts.items():
            if key not in DEFAULT_FIT_OPTS:
                logging.warning('Ignoring unknown fit option %r', key)
            else:
                merged[key] = value
    return merged


def robust_scale(res):
    '''
    Robust residual scale: 1.4826 * median absolute deviation (MAD).

    Returns 0.0 when MAD is zero (e.g. too few points or all-equal
    residuals), so the caller can detect the degenerate case and fall back.

    Input:
        res = 1-D array of residuals

    Output:
        scalar robust standard-deviation estimate
    '''
    res = np.asarray(res, dtype=float)
    if res.size == 0:
        return 0.0
    med = np.median(res)
    return 1.4826 * np.median(np.abs(res - med))


def flag_outliers(res, sigma, edt_mask, nsigma, scale_mode='mad',
                  nrms=1.0, restore_factor=0.9):
    '''
    One n-sigma editing pass, mirroring tsfit's edit_ns but with a
    configurable scale.

    The editor's own flags live in ``edt_mask``. A point is flagged when it is
    not currently editor-flagged and its normalised residual exceeds nsigma;
    an already-flagged point is restored when its normalised residual drops to
    <= restore_factor * nsigma (hysteresis, to avoid flip-flopping near the
    threshold). Points excluded by other, earlier masks (e.g. max_sigma) are
    not touched here -- the caller keeps those masks separate.

    Input:
        res            = 1-D residuals (same length as sigma)
        sigma          = 1-D formal sigmas
        edt_mask       = bool mask of the editor's current flags
        nsigma         = rejection threshold (units of the chosen scale)
        scale_mode     = 'mad' (robust residual scale) or 'sigma' (sigma*nrms)
        nrms           = white-noise NRMS, only used when scale_mode == 'sigma'
        restore_factor = hysteresis multiplier for restoring points

    Output:
        new_edt_mask = bool mask with the same length as ``edt_mask``
    '''
    edt_mask = np.asarray(edt_mask, dtype=bool)
    res = np.asarray(res, dtype=float)
    sigma = np.asarray(sigma, dtype=float)

    new_edt_mask = edt_mask.copy()

    if scale_mode == 'mad':
        scale = robust_scale(res[~edt_mask])
        if scale <= 0.0:
            logging.warning('robust_scale is zero; skipping this editing pass')
            return new_edt_mask
        denom = np.full(res.shape, scale)
    elif scale_mode == 'sigma':
        denom = sigma * nrms
    else:
        logging.warning('Unknown outlier_scale %r; skipping editing', scale_mode)
        return new_edt_mask

    # Points with non-positive sigma cannot be normalised by the sigma scale;
    # treat them as infinite residual so they are (conservatively) flagged.
    safe_denom = np.where(denom > 0.0, denom, 1.0)
    err = np.abs(res / safe_denom)
    err[denom <= 0.0] = np.inf

    # Flag currently-unflagged points beyond the threshold.
    flag = (~edt_mask) & (err > nsigma)
    # Restore previously-flagged points now back within the hysteresis band.
    restore = edt_mask & (err <= nsigma * restore_factor)

    new_edt_mask[flag] = True
    new_edt_mask[restore] = False
    return new_edt_mask


def realistic_sigma(t, res, sigma, good=None):
    '''
    Herring bin-averaging + first-order Gauss-Markov extrapolation, ported
    from GAMIT/GLOBK gen_util/real_stats.f.

    The returned scale is unit-independent (it is a dimensionless ratio of
    empirical bin-mean variance to the white-noise expectation), so both
    millimetre and metre inputs give the same result.

    Input:
        t     = times in days (only differences matter; not required to be
                sorted, but sorted internally)
        res   = post-fit residuals
        sigma = formal sigmas
        good  = optional bool mask of points to use (all True when None)

    Output:
        (sig_scale, taufin) with taufin in days, or (None, None) when the
        series is too short / degenerate for the FOGM fit
    '''
    t = np.asarray(t, dtype=float)
    res = np.asarray(res, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    if good is None:
        good = np.ones(t.shape, dtype=bool)
    else:
        good = np.asarray(good, dtype=bool)

    # Use only finite, positive-sigma, good points.
    use = good & np.isfinite(t) & np.isfinite(res) & (sigma > 0.0)
    t = t[use]
    res = res[use]
    sigma = sigma[use]

    if t.size < 2:
        return None, None

    # real_stats assumes sorted times; normalise the origin so the first
    # epoch is zero (differences only).
    order = np.argsort(t)
    t = t[order]
    res = res[order]
    sigma = sigma[order]
    t = t - t[0]

    stt = t[0]
    ent = t[-1]
    span = ent - stt

    minav = 7.0
    maxav = span / 10.0
    if maxav < minav * 3.0:
        minav = int(maxav / 4.0) + 1
    if minav < 1.0:
        return None, None

    numav = int(maxav / minav)
    if numav > _MAX_AV:
        logging.warning('realistic_sigma: numav %d clipped to %d', numav, _MAX_AV)
        numav = _MAX_AV
    if numav < 2:
        return None, None

    chi2 = np.empty(numav, dtype=float)
    tims = np.empty(numav, dtype=float)
    inv_sigma2 = 1.0 / sigma ** 2

    for i in range(numav):
        dt = (i + 1) * minav
        sumavs = 0.0
        num = 0
        nwin = int(np.rint((ent - stt) / dt)) + 1
        for it in range(nwin):
            lo = stt + it * dt
            hi = lo + dt
            in_bin = (t >= lo) & (t < hi)
            if not in_bin.any():
                continue
            w = inv_sigma2[in_bin]
            sw = w.sum()
            if sw <= 0.0:
                continue
            wmean = np.sum(res[in_bin] * w) / sw
            varm = 1.0 / sw
            sumavs += wmean * wmean / varm
            num += 1
        if num == 0:
            # No usable bin at this averaging interval: cannot form a chi2
            # curve, so the FOGM fit is not attempted (real_stats.f divides
            # by num here and produces NaN; we fail cleanly instead).
            return None, None
        chi2[i] = sumavs / num
        tims[i] = dt

    # Scan the trial correlation times and keep the best FOGM curve.
    rmmin = np.inf
    sig_scale = None
    taufin = None
    for tau in _TAUS:
        ef = 1.0 - np.exp(-tims / tau)
        # real_stats computes alpha as the mean of chi2/ef (ratio estimator),
        # not a least-squares fit. Guard against ef == 0 for tiny tims/tau.
        nz = ef > 0.0
        if not nz.all():
            continue
        alpha = np.mean(chi2 / ef)
        rmsum = np.sum((chi2 - ef * alpha) ** 2)
        if rmsum < rmmin:
            rmmin = rmsum
            sig_scale = alpha
            taufin = tau

    if sig_scale is None or not np.isfinite(sig_scale) or sig_scale < 0.0:
        return None, None
    return float(np.sqrt(sig_scale)), float(taufin)
