'''
Tests for realistic_sigma: the Herring bin-averaging + FOGM port from
GAMIT/GLOBK gen_util/real_stats.f.

The first test pins the exact (sig_scale, taufin) values returned by the
Fortran routine for simple deterministic inputs, cross-checked against a
gfortran-compiled copy of the original real_stats.f.
'''
import numpy as np
import pytest

from pytsfit.qualitycontrol import realistic_sigma, merge_fit_opts


def _series(npts, step_days):
    '''Simple deterministic residual: sin + a gentle ramp (matching the
    Fortran cross-check driver).'''
    t = np.arange(npts) * step_days
    res = np.sin(0.03 * t) + 0.01 * t
    sigma = np.ones(npts)
    return t, res, sigma


def test_matches_fortran_reference():
    # (npts, step_days, sig_scale, taufin) — values produced by real_stats.f
    # compiled from /home/zhao/gg/kf/gen_util/real_stats.f.
    cases = [
        (100,  1, 7.208674523500, 32.0),
        (300,  1, 45.735318989633, 512.0),
        (50,   7, 13.295705813192, 256.0),
        (500,  3, 160.779263845764, 1024.0),
        (1200, 3, 765.743105658203, 4096.0),
        (80,   7, 20.478003313000, 256.0),
    ]
    for npts, step, want_scale, want_tau in cases:
        t, res, sigma = _series(npts, step)
        got_scale, got_tau = realistic_sigma(t, res, sigma)
        assert got_tau == want_tau, (npts, step, got_tau)
        assert got_scale == pytest.approx(want_scale, rel=1e-12, abs=1e-12)


def test_colored_noise_has_larger_scale_than_white_noise():
    rng = np.random.default_rng(20240828)
    npts, step = 800, 3
    t = np.arange(npts) * step
    sigma = np.ones(npts)

    white = rng.normal(0.0, 1.0, npts)
    rw = np.cumsum(rng.normal(0.0, 0.3, npts))    # random-walk / FOGM-like

    ws, _ = realistic_sigma(t, white, sigma)
    cs, ctau = realistic_sigma(t, rw, sigma)
    # Coloured noise needs a larger scale than white noise, and a non-trivial
    # correlation time.
    assert cs > ws
    assert ctau >= 32.0


def test_too_short_series_returns_none():
    # Two epochs is too short to form >= 2 averaging intervals.
    t, res, sigma = _series(2, 7)
    assert realistic_sigma(t, res, sigma) == (None, None)
    # Three epochs is also too short.
    t, res, sigma = _series(3, 7)
    assert realistic_sigma(t, res, sigma) == (None, None)


def test_unsorted_input_matches_sorted():
    t, res, sigma = _series(500, 3)
    rng = np.random.default_rng(1)
    perm = rng.permutation(t.size)
    ss_sorted, tau_sorted = realistic_sigma(t, res, sigma)
    ss_perm, tau_perm = realistic_sigma(t[perm], res[perm], sigma[perm])
    assert tau_perm == tau_sorted
    assert ss_perm == pytest.approx(ss_sorted, rel=1e-12)


def test_scale_is_dimensionless():
    # The result is a dimensionless ratio, so scaling both residual and sigma
    # (e.g. metres -> millimetres) must not change it.
    t, res, sigma = _series(600, 3)
    ss1, tau1 = realistic_sigma(t, res, sigma)
    ss2, tau2 = realistic_sigma(t, res * 1000.0, sigma * 1000.0)
    assert tau1 == tau2
    assert ss1 == pytest.approx(ss2, rel=1e-12)


def test_merge_fit_opts_defaults_and_override():
    defaults = merge_fit_opts(None)
    assert defaults['sigma_scale'] == 'nrms'
    assert defaults['outlier'] is False
    assert defaults['min_rsig'] == 30

    merged = merge_fit_opts({'sigma_scale': 'realistic', 'unknown_key': 1})
    assert merged['sigma_scale'] == 'realistic'
    assert 'unknown_key' not in merged
    # merging must not mutate the shared defaults dict
    assert merge_fit_opts(None)['sigma_scale'] == 'nrms'
