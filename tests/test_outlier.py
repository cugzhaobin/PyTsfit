'''
Tests for the iterative outlier editor and its wiring into tsfitting.

The editor mirrors tsfit's edit_ns: fit -> flag by n-sigma -> refit, until the
flagged set stops changing, with a restore_factor hysteresis band.
'''
import numpy as np
import pytest

from pytsfit import tsfitting, correction, flag_outliers, robust_scale, merge_fit_opts


def _series(vel_true=5.0, noise=1.0, npts=500, seed=7):
    '''Clean constant + linear + annual series with known truth.'''
    rng = np.random.default_rng(seed)
    t = 2010.0 + np.arange(npts) / 365.25
    obs = (10.0
           + vel_true * (t - t.mean())
           + 2.0 * np.sin(2 * np.pi * t) + 1.0 * np.cos(2 * np.pi * t)
           + rng.normal(0.0, noise, npts))
    sigma = np.full(npts, noise)
    return t, obs, sigma


def _param_dict():
    return {'constant': True, 'linear': True, 'eqlist': [], 'eqpostlist': [],
            'brklist': [], 'ANN': True, 'SANN': False,
            'correct': correction('', '', '')}


def _run(t, obs, sigma, fit_opts):
    pd = _param_dict()
    run = tsfitting('TEST', 103.2, 35.5, t, obs, sigma, pd, 'N', fit_opts=fit_opts)
    run.doFitting()
    return run


# ---------------------------------------------------------------------------
# flag_outliers / robust_scale: direct unit tests
# ---------------------------------------------------------------------------
def test_flag_outliers_flags_and_restores_with_hysteresis():
    res = np.array([0.1, 0.2, 20.0, 0.3, 0.4])      # index 2 is a spike
    sigma = np.ones(5)
    edt = np.zeros(5, dtype=bool)
    edt = flag_outliers(res, sigma, edt, nsigma=4.0, scale_mode='mad')
    assert edt.tolist() == [False, False, True, False, False]

    # Hysteresis is cleanest to pin in 'sigma' mode, where the threshold is in
    # units of residual/sigma*nrms and so directly comparable to nsigma.
    res = np.array([0.1, 0.2, 20.0, 0.3, 0.4])
    edt = flag_outliers(res, sigma, np.zeros(5, dtype=bool), nsigma=4.0,
                        scale_mode='sigma', nrms=1.0)
    assert edt.tolist() == [False, False, True, False, False]

    res2 = res.copy()
    res2[2] = 3.5                                  # <= 0.9 * 4.0 -> restore
    edt2 = flag_outliers(res2, sigma, edt, nsigma=4.0, scale_mode='sigma',
                         nrms=1.0, restore_factor=0.9)
    assert edt2.tolist() == [False, False, False, False, False]

    res3 = res.copy()
    res3[2] = 3.7                                  # > 0.9 * 4.0 -> keep flagged
    edt3 = flag_outliers(res3, sigma, edt, nsigma=4.0, scale_mode='sigma',
                         nrms=1.0, restore_factor=0.9)
    assert edt3.tolist() == [False, False, True, False, False]


def test_robust_scale_matches_expected_mad():
    res = np.array([1.0, 2.0, 3.0, 4.0, 100.0])
    # median = 3, abs deviations median = 1 -> 1.4826
    assert robust_scale(res) == pytest.approx(1.4826)


def test_merge_fit_opts_returns_defaults():
    assert merge_fit_opts(None)['outlier'] is False
    assert merge_fit_opts({'outlier': True, 'nsigma': 6.0})['nsigma'] == 6.0


# ---------------------------------------------------------------------------
# End-to-end editing through tsfitting
# ---------------------------------------------------------------------------
def test_outlier_editor_removes_injected_spikes():
    t, obs, sigma = _series()
    npts = obs.size
    outlier_idx = np.array([50, 200, 400])
    obs[outlier_idx] += 40.0                       # ~40 sigma spikes

    run = _run(t, obs, sigma, {'outlier': True, 'nsigma': 4.0, 'outlier_scale': 'mad'})
    assert run.good.dtype == bool
    # exactly the injected positions are editor-flagged
    assert set(np.where(run.edt_mask)[0].tolist()) == set(outlier_idx.tolist())
    assert run.nedit == outlier_idx.size


def test_outlier_editor_recovers_true_velocity():
    t, obs, sigma = _series(vel_true=5.0)
    obs[[80, 300, 450]] += 50.0
    run = _run(t, obs, sigma, {'outlier': True, 'nsigma': 4.0, 'outlier_scale': 'mad'})
    vel_idx = np.where(run.flag2 == 'VELOCITY')[0][0]
    assert run.param[vel_idx] == pytest.approx(5.0, abs=0.15)
    assert run.niter <= 10


def test_components_are_independent():
    t, obs, sigma = _series()
    bad = 40.0
    u_obs = obs.copy()
    u_obs[100] += bad

    # North (clean) and Up (with a spike) are fitted independently.
    nrun = _run(t, obs, sigma, {'outlier': True, 'nsigma': 4.0, 'outlier_scale': 'mad'})
    urun = _run(t, u_obs, sigma, {'outlier': True, 'nsigma': 4.0, 'outlier_scale': 'mad'})
    assert nrun.nedit == 0
    assert nrun.good.all()
    assert urun.nedit == 1
    assert not urun.good[100]


def test_outlier_disabled_leaves_everything_clean_and_unchanged():
    t, obs, sigma = _series()
    obs[[10, 100]] += 40.0
    run = _run(t, obs, sigma, {'outlier': False})
    assert run.good.all()
    assert run.nedit == 0


def test_sigma_scale_mode_uses_formal_sigma():
    t, obs, sigma = _series(noise=1.0)
    obs[300] += 30.0
    # sigma*nrms ~ 1; the spike is ~30 sigma -> flagged with mode='sigma'.
    run = _run(t, obs, sigma, {'outlier': True, 'nsigma': 4.0, 'outlier_scale': 'sigma'})
    assert run.edt_mask[300]


def test_max_sigma_prefilter_flags_large_sigma_points():
    t, obs, sigma = _series()
    sigma = sigma.copy()
    sigma[[5, 6, 7]] = 100.0
    run = _run(t, obs, sigma, {'outlier': False, 'max_sigma': 20.0})
    assert set(np.where(~run.good)[0].tolist()) == {5, 6, 7}
