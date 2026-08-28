'''
Regression tests for the specific bugs fixed alongside the module split.

These pin the *fixed* behaviour. Several of them would have failed before the
fix commit (e.g. constructing two catalogs returned a duplicated list, and the
BREAK branch of setBoundAndInit raised IndexError on offsetdata[:,4]).
'''
import numpy as np
import pytest

from pytsfit.PyTsfit import (neuData, posData, earthquake, eqcatalog, eqPost,
                             eqPostList, offset, breakcatalog, correction,
                             tsfitting, output_postseismic_disp)

import fixtures


# ---------------------------------------------------------------------------
# Catalog instance independence
# ---------------------------------------------------------------------------
def test_eqcatalog_instances_do_not_share_list(tmp_path):
    p = fixtures.write_eq_rename(str(tmp_path / 'eq_rename'))
    a = eqcatalog(p)
    b = eqcatalog(p)
    assert len(a.eqlist) == 3
    assert len(b.eqlist) == 3
    assert a.eqlist is not b.eqlist
    # mutating one must not leak into the other
    a.eqlist.pop()
    assert len(a.eqlist) == 2
    assert len(b.eqlist) == 3


def test_breakcatalog_instances_do_not_share_list(tmp_path):
    p = fixtures.write_eq_rename(str(tmp_path / 'eq_rename'))
    a = breakcatalog(p)
    b = breakcatalog(p)
    assert len(a.breaklist) == 3
    assert len(b.breaklist) == 3
    assert a.breaklist is not b.breaklist


def test_eqpostlist_instances_do_not_share_list(tmp_path):
    p = fixtures.write_eq_rename(str(tmp_path / 'eq_rename'))
    eq = eqcatalog(p)
    a = eqPostList(p, eq)
    b = eqPostList(p, eq)
    assert len(a.eqpostlist) == 1
    assert len(b.eqpostlist) == 1
    assert a.eqpostlist is not b.eqpostlist


def test_eq_exp_parsed_as_exp(tmp_path):
    txt = ('  eq_def TA 35.60 103.30 800 12 2010 4 13 12 0\n'
           '  eq_exp TA 20 50\n')
    p = tmp_path / 'eq_rename'
    p.write_text(txt)
    eq = eqcatalog(str(p))
    post = eqPostList(str(p), eq)
    assert post.eqpostlist[0].method == 'EXP'


# ---------------------------------------------------------------------------
# offset.fun_offset: self.decyr is a float, not a callable
# ---------------------------------------------------------------------------
def test_offset_fun_offset():
    off = offset([2010, 4, 13, 12, 0], 'TEST')
    t = np.array([2009.0, 2010.3, 2011.0])
    y = off.fun_offset(t, 5.0)
    assert np.allclose(y, [0.0, 5.0, 5.0])


# ---------------------------------------------------------------------------
# setBoundAndInit BREAK branch: read the year from offsetyear, not offsetdata[:,4]
# ---------------------------------------------------------------------------
def test_break_bounds_use_offsetyear(tmp_path):
    velfile, offsetfile, periodfile = fixtures.write_correction(tmp_path)
    cor = correction(velfile, offsetfile, periodfile)

    eqfile = tmp_path / 'eq_rename'
    fixtures.write_eq_rename(str(eqfile))
    bk = breakcatalog(str(eqfile))

    t = 2010.0 + np.arange(1200) / 365.25
    t0 = bk.breaklist[0].decyr
    obs = 10.0 + 2.0 * (t - t.mean()) + 3.0 * np.heaviside(t - t0, 0.0)
    sig = np.full(t.size, 1.0)

    param_dict = {'constant': True, 'linear': True, 'eqlist': [], 'eqpostlist': [],
                  'brklist': bk.breaklist, 'ANN': False, 'SANN': False,
                  'correct': cor}
    run = tsfitting('TEST', fixtures.LON, fixtures.LAT, t, obs, sig, param_dict, 'N')

    lb, ub, pinit = run.setBoundAndInit()   # must not raise IndexError
    assert len(lb) == run.nparam
    assert np.isfinite(pinit).all()


# ---------------------------------------------------------------------------
# Mutable default arguments
# ---------------------------------------------------------------------------
def _make_run():
    t = np.linspace(2010.0, 2019.0, 500)
    obs = 10.0 + 2.0 * (t - t.mean())
    sig = np.full(t.size, 1.0)
    param_dict = {'constant': True, 'linear': True, 'eqlist': [], 'eqpostlist': [],
                  'brklist': [], 'ANN': False, 'SANN': False,
                  'correct': correction('', '', '')}
    run = tsfitting('TEST', 103.2, 35.5, t, obs, sig, param_dict, 'N')
    run.doFitting()
    return run


def test_get_mod_get_correct_defaults_are_not_shared():
    r1 = _make_run()
    mt1, m1 = r1.get_mod()
    r2 = _make_run()
    mt2, m2 = r2.get_mod()
    assert mt1.size == mt2.size
    # full-range default is independent of any prior time_span call
    r1.get_mod([2012.0, 2015.0])
    mt3, _ = r1.get_mod()
    assert mt3.size == mt1.size


def test_tsfitting_time_range_default_is_not_shared():
    t = np.linspace(2010.0, 2019.0, 500)
    obs = 10.0 + 2.0 * (t - t.mean())
    sig = np.full(t.size, 1.0)
    pd = {'constant': True, 'linear': True, 'eqlist': [], 'eqpostlist': [],
          'brklist': [], 'ANN': False, 'SANN': False, 'correct': correction('', '', '')}
    r1 = tsfitting('TEST', 103.2, 35.5, t, obs, sig, pd, 'N')
    r2 = tsfitting('TEST', 103.2, 35.5, t, obs, sig, pd, 'N', [2013.0, 2016.0])
    assert r1.t.size == t.size
    assert r2.t.size < t.size


# ---------------------------------------------------------------------------
# tsfitting.plot_obs_mod: self.parm -> self.param
# ---------------------------------------------------------------------------
def test_tsfitting_plot_obs_mod_uses_param(monkeypatch):
    import matplotlib.pyplot as plt
    run = _make_run()
    calls = []
    # Stub the plotting calls: the bug we are pinning is the self.parm ->
    # self.param fix; the (separate) model-vs-time length mismatch is not this
    # test's concern.
    monkeypatch.setattr(plt, 'plot', lambda *a, **k: calls.append(a))
    monkeypatch.setattr(plt, 'title', lambda *a, **k: None)
    monkeypatch.setattr(plt, 'xlabel', lambda *a, **k: None)
    monkeypatch.setattr(plt, 'ylabel', lambda *a, **k: None)
    monkeypatch.setattr(plt, 'show', lambda *a, **k: None)
    # would raise AttributeError on the old self.parm reference
    run.plot_obs_mod()
    assert calls, 'plot was never reached'
    # the model was built from run.param (a float array), not a stale attribute
    assert np.asarray(run.param).size == run.nparam


# ---------------------------------------------------------------------------
# output_postseismic_disp: format-string index and unreachable return
# ---------------------------------------------------------------------------
def test_output_postseismic_disp_no_index_error():
    # A run whose time series predates the requested time_span. The empty-model
    # guard then runs; with the old "{1:s}".format(site) it would raise
    # IndexError on the warning path.
    t = np.linspace(2008.0, 2010.0, 300)
    obs = 10.0 + 2.0 * (t - t.mean())
    sig = np.full(t.size, 1.0)
    param_dict = {'constant': True, 'linear': True, 'eqlist': [], 'eqpostlist': [],
                  'brklist': [], 'ANN': False, 'SANN': False,
                  'correct': correction('', '', '')}
    run = tsfitting('TEST', 103.2, 35.5, t, obs, sig, param_dict, 'N')
    run.doFitting()
    # time_span fully after the data: get_correct yields no observed points
    output_postseismic_disp(run, run, run, [2015.0, 2016.0], fid=None)
