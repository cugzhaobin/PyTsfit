'''
Build a complete, JSON-serializable snapshot of PyTsfit's observable behaviour.

The snapshot is the contract the refactor must preserve. Generate it once
against the pre-refactor code (tests/golden/baseline.json), then re-generate it
after moving code around; any difference is a behaviour change.

Everything is imported through ``pytsfit.PyTsfit`` on purpose: that name is the
public surface today and stays a compatibility shim after the split, so the
same snapshot code runs unchanged on both sides of the refactor.
'''
import json
import os

import matplotlib
matplotlib.use('Agg')          # no display in CI; must precede pyplot import
import numpy as np

from pytsfit import PyTsfit as P

import fixtures


def _f(x):
    '''Coerce a numpy scalar to a plain float that json round-trips exactly.'''
    return float(x)


def _arr(a, limit=None):
    '''Serialize an array as a list of exact floats, optionally truncated.'''
    a = np.asarray(a, dtype=float).ravel()
    if limit is not None:
        a = a[:limit]
    return [_f(v) for v in a]


def _stats(a):
    '''Compact fingerprint of a large array: length plus exact moments.'''
    a = np.asarray(a, dtype=float).ravel()
    if a.size == 0:
        return {'n': 0}
    return {'n': int(a.size), 'sum': _f(a.sum()), 'min': _f(a.min()),
            'max': _f(a.max()), 'first5': _arr(a, 5), 'last5': _arr(a[-5:])}


def _reset_catalog_class_state():
    '''
    eqcatalog/breakcatalog/eqPostList currently accumulate into class-level
    lists, so results depend on how many instances were built earlier in the
    process. Clear them so the snapshot is reproducible regardless.
    Harmless once those become per-instance attributes.
    '''
    for cls, attr in ((P.eqcatalog, 'eqlist'),
                      (P.breakcatalog, 'breaklist'),
                      (P.eqPostList, 'eqpostlist')):
        if isinstance(getattr(cls, attr, None), list):
            setattr(cls, attr, [])


def _snapshot_time_and_geo():
    '''Lock down the GPSTime / geotools functions PyTsfit depends on.'''
    gpstime = P.gpstime
    gt = P.gt

    mjds = list(np.arange(50000.0, 61000.0, 250.0)) + [55197.25, 55299.5, 58000.125]
    decyrs = [_f(gpstime.jd_to_decyrs(m)) for m in mjds]

    dates = [[1999, 1, 1, 0, 0], [2008, 5, 12, 6, 28], [2011, 3, 11, 5, 46],
             [2016, 2, 29, 23, 59], [2020, 12, 31, 12, 0], [2024, 6, 15, 7, 31]]
    jds = [_f(gpstime.ymdhms_to_jd(list(d), 0)) for d in dates]

    rng = np.random.default_rng(7)
    xy = []
    for _ in range(50):
        lat, lon = rng.uniform(18, 53), rng.uniform(73, 135)
        olat, olon = rng.uniform(18, 53), rng.uniform(73, 135)
        v = gt.llh2localxy([lat, lon], [olat, olon, 10.0])
        xy.append([_f(v[0]), _f(v[1])])

    return {'jd_to_decyrs': decyrs, 'ymdhms_to_jd': jds, 'llh2localxy': xy}


def _snapshot_readers(tmpdir):
    '''posData and neuData parsing.'''
    out = {}

    posfile = str(tmpdir / 'TEST.pos')
    fixtures.write_pos(posfile)
    pos = P.posData(posfile)
    out['posData'] = {
        'site': pos.site, 'lat': _f(pos.lat), 'lon': _f(pos.lon), 'hei': _f(pos.hei),
        'decyr': _stats(pos.decyr), 'N': _stats(pos.N), 'E': _stats(pos.E),
        'U': _stats(pos.U), 'SN': _stats(pos.SN), 'SE': _stats(pos.SE), 'SU': _stats(pos.SU),
    }

    neufile = str(tmpdir / 'TEST.neu')
    fixtures.write_neu(neufile)
    neu = P.neuData(neufile)
    out['neuData'] = {
        'site': neu.site, 'lat': _f(neu.lat), 'lon': _f(neu.lon),
        'decyr': _stats(neu.decyr), 'N': _stats(neu.N), 'E': _stats(neu.E), 'U': _stats(neu.U),
    }
    return out


def _snapshot_catalogs(eqfile):
    '''eqcatalog / breakcatalog / eqPostList parsing and selection.'''
    _reset_catalog_class_state()
    eq = P.eqcatalog(eqfile)
    bk = P.breakcatalog(eqfile)
    ep = P.eqPostList(eqfile, eq)

    return {
        'eqcatalog': [{'code': e.code, 'location': _arr(e.location),
                       'epoch': [int(v) for v in e.epoch],
                       'distance': _f(e.distance), 'decyr': _f(e.decyr)}
                      for e in eq.eqlist],
        'getEQ_TA_decyr': _f(eq.getEQ('TA').decyr),
        'breakcatalog': [{'site': b.site, 'epoch': [int(v) for v in b.epoch],
                          'decyr': _f(b.decyr)} for b in bk.breaklist],
        'eqPostList': [{'code': q.eq.code, 'method': q.method,
                        'mintau': str(q.mintau), 'maxtau': str(q.maxtau)}
                       for q in ep.eqpostlist],
    }, eq, bk, ep


def _snapshot_correction(tmpdir):
    '''class correction, with files present and with files absent.'''
    velfile, offsetfile, periodfile = fixtures.write_correction(tmpdir)

    present = P.correction(velfile, offsetfile, periodfile)
    absent = P.correction('', '', '')
    no_offset = P.correction(velfile, '', periodfile)

    return {
        'present': {
            'correct_velo': bool(present.correct_velo),
            'correct_offset': bool(present.correct_offset),
            'correct_period': bool(present.correct_period),
            'velsite': [str(s) for s in np.atleast_1d(present.velsite)],
            'veldata': _arr(present.veldata),
            'offsetsite': [str(s) for s in np.atleast_1d(present.offsetsite)],
            'offsetdata': _arr(present.offsetdata),
            'offsetyear': _arr(present.offsetyear),
            'periodsite': [str(s) for s in np.atleast_1d(present.periodsite)],
            'perioddata': _arr(present.perioddata),
        },
        'absent': {
            'correct_velo': bool(absent.correct_velo),
            'correct_offset': bool(absent.correct_offset),
            'correct_period': bool(absent.correct_period),
        },
    }, present, absent, no_offset


def _fit_case(pos, param_dict, component, time_range):
    '''Run one tsfitting configuration and capture everything it exposes.'''
    obs = {'N': (pos.N, pos.SN), 'E': (pos.E, pos.SE), 'U': (pos.U, pos.SU)}[component]
    run = P.tsfitting(pos.site, pos.lon, pos.lat, pos.decyr, obs[0], obs[1],
                      param_dict, component, time_range)

    rec = {
        'nparam': int(run.nparam),
        'flag': [str(v) for v in run.flag],
        'flag2': [str(v) for v in run.flag2],
        'ieqlist': [e.code for e in run.ieqlist],
        'ieqpostlist': [q.eq.code for q in run.ieqpostlist],
        'ibrklist': [_f(b.decyr) for b in run.ibrklist],
        'midt': _f(run.midt) if run.midt is not None else None,
        'n_obs': int(np.size(run.t)),
    }

    if run.nparam:
        lb, ub, pinit = run.setBoundAndInit()
        rec['bounds'] = {'lb': [_f(v) for v in lb], 'ub': [_f(v) for v in ub],
                         'pinit': [_f(v) for v in pinit]}

    popt = run.doFitting()
    rec['popt'] = _arr(popt)
    if popt.size:
        rec['cov_diag'] = _arr(np.diag(run.cov))
        rec['wrms'] = _f(run.wrms)
        rec['res'] = _stats(run.res)

        mt, m = run.get_mod()
        rec['get_mod'] = {'t': _stats(mt), 'm': _stats(m)}
        mt2, m2 = run.get_mod([2012.0, 2016.0])
        rec['get_mod_span'] = {'t': _stats(mt2), 'm': _stats(m2)}

        mod_dict = {'detrend': True, 'debreak': True, 'deeqoffset': False,
                    'depost': False, 'deseason': False}
        oc, mc = run.get_correct(mod_dict)
        rec['get_correct'] = {'obs': _stats(oc), 'mod': _stats(mc)}

        mod_dict2 = {'detrend': False, 'debreak': False, 'deeqoffset': False,
                     'depost': True, 'deseason': False}
        oc2, mc2 = run.get_correct(mod_dict2, time_span=[2012.0, 2016.0], eqcode='TA')
        rec['get_correct_span'] = {'obs': _stats(oc2), 'mod': _stats(mc2)}

    return rec, run


def _snapshot_fits(pos, eq, bk, ep, cor_present, cor_absent, cor_no_offset):
    '''Several tsfitting configurations, each in all three components.'''
    base_unconstrained = {
        'constant': True, 'linear': True,
        'eqlist': eq.eqlist, 'eqpostlist': [], 'brklist': bk.breaklist,
        'ANN': True, 'SANN': True, 'correct': cor_absent,
    }
    cases = {
        'full_unconstrained': (dict(base_unconstrained), [-np.inf, np.inf]),
        'with_postseismic': (dict(base_unconstrained, eqpostlist=ep.eqpostlist), [-np.inf, np.inf]),
        'constrained': (dict(base_unconstrained, correct=cor_no_offset), [-np.inf, np.inf]),
        'linear_only': (dict(base_unconstrained, eqlist=[], brklist=[],
                             ANN=False, SANN=False), [-np.inf, np.inf]),
        'time_window': (dict(base_unconstrained), [2012.0, 2017.0]),
    }

    out = {}
    runs = {}
    for name, (param_dict, time_range) in cases.items():
        out[name] = {}
        runs[name] = {}
        for comp in ('N', 'E', 'U'):
            rec, run = _fit_case(pos, param_dict, comp, time_range)
            out[name][comp] = rec
            runs[name][comp] = run
    return out, runs


def _capture(fn, *args, **kwargs):
    '''Run an output_* writer against a file handle and return the text.'''
    path = 'capture.tmp'
    with open(path, 'w') as fid:
        fn(*args, fid=fid, **kwargs)
    with open(path) as fid:
        text = fid.read()
    os.remove(path)
    return text


def _snapshot_outputs(runs, tmpdir):
    '''Exact text emitted by every output_* function.'''
    nrun = runs['with_postseismic']['N']
    erun = runs['with_postseismic']['E']
    urun = runs['with_postseismic']['U']

    out = {}
    for fmt in ('GMT', 'IOS3D', 'GLOBK', 'DETAIL'):
        out['velo_' + fmt] = _capture(P.output_velo, nrun, erun, urun, fmt=fmt)
    for fmt in ('GMT2D', 'GMT3D', 'IOS3D'):
        out['eqoffset_' + fmt] = _capture(P.output_eqoffset, nrun, erun, urun, fmt=fmt)

    out['postseismic_velo'] = _capture(P.output_postseismic_velo, nrun, erun, urun,
                                       [2011.0, 2016.0])
    out['postseismic_disp'] = _capture(P.output_postseismic_disp, nrun, erun, urun,
                                       [2011.0, 2016.0])
    out['summary'] = _capture(P.output_summary, nrun, erun, urun)

    # output_period appends to a fixed filename in the working directory.
    period_path = 'period.dat'
    if os.path.exists(period_path):
        os.remove(period_path)
    P.output_period(nrun, erun, urun, nrun.param, erun.param, urun.param)
    with open(period_path) as fid:
        out['period'] = fid.read()
    os.remove(period_path)

    # output_obs_mod writes <site>_obs.dat and <site>_mod.dat.
    mod_dict = {'detrend': True, 'debreak': True, 'deeqoffset': False,
                'depost': False, 'deseason': False}
    P.output_obs_mod(nrun, erun, urun, nrun.param, erun.param, urun.param, mod_dict)
    for suffix in ('obs', 'mod'):
        name = '{}_{}.dat'.format(nrun.site, suffix)
        with open(name) as fid:
            text = fid.read()
        # Files are large; fingerprint the header plus first and last records.
        lines = text.splitlines()
        out['obs_mod_' + suffix] = {'n_lines': len(lines), 'head': lines[:4],
                                    'tail': lines[-2:]}
        os.remove(name)

    # output_postseismic_ts writes <site>_<eqcode>.neu.
    P.output_postseismic_ts(nrun, erun, urun, [2011.0, 2016.0], eqcode='TA')
    ts_name = '{}_TA.neu'.format(nrun.site)
    with open(ts_name) as fid:
        lines = fid.read().splitlines()
    out['postseismic_ts'] = {'n_lines': len(lines), 'head': lines[:4], 'tail': lines[-2:]}
    os.remove(ts_name)

    # output_param writes <site>_par.npy.
    P.output_param(nrun, erun, urun, nrun.param, erun.param, urun.param)
    par = np.load('{}_par.npy'.format(nrun.site), allow_pickle=True).item()
    out['param'] = {
        'site': par['site'],
        'Vn': _arr(par['Vn']), 'Sn': _arr(par['Sn']),
        'Ve': _arr(par['Ve']), 'Se': _arr(par['Se']),
        'Vu': _arr(par['Vu']), 'Su': _arr(par['Su']),
        'eqcode': list(par['eqinfo']['eqcode']) if par['eqinfo'] else None,
        'eqyear': _arr(par['eqinfo']['eqyear']) if par['eqinfo'] else None,
        'eqdisp': [_arr(d) for d in par['eqinfo']['eqdisp']] if par['eqinfo'] else None,
    }
    os.remove('{}_par.npy'.format(nrun.site))

    # plot_obs_mod is exercised for crash-freedom only; pixels are not a contract.
    plot_dict = {'detrend': True, 'debreak': True, 'deeqoffset': False,
                 'depost': False, 'deseason': False,
                 'figformat': 'jpg', 'showfig': False}
    P.plot_obs_mod(nrun, erun, urun, nrun.param, erun.param, urun.param, plot_dict)
    fig_name = '{}.jpg'.format(nrun.site)
    out['plot_obs_mod_ran'] = os.path.exists(fig_name)
    if os.path.exists(fig_name):
        os.remove(fig_name)

    return out


def build_snapshot(tmpdir):
    '''
    Build the full snapshot. ``tmpdir`` is a pathlib.Path that this function
    also chdir's into, because several output_* writers use fixed filenames
    relative to the working directory.
    '''
    cwd = os.getcwd()
    os.chdir(str(tmpdir))
    try:
        snap = {}
        snap['time_and_geo'] = _snapshot_time_and_geo()
        snap['readers'] = _snapshot_readers(tmpdir)

        eqfile = str(tmpdir / 'eq_rename')
        fixtures.write_eq_rename(eqfile)
        snap['catalogs'], eq, bk, ep = _snapshot_catalogs(eqfile)

        snap['correction'], cor_present, cor_absent, cor_no_offset = _snapshot_correction(tmpdir)

        pos = P.posData(str(tmpdir / 'TEST.pos'))
        snap['fits'], runs = _snapshot_fits(pos, eq, bk, ep, cor_present, cor_absent, cor_no_offset)
        snap['outputs'] = _snapshot_outputs(runs, tmpdir)
        return snap
    finally:
        os.chdir(cwd)


def main():
    import argparse
    import pathlib
    import tempfile

    parser = argparse.ArgumentParser(description='Write a PyTsfit behaviour snapshot.')
    parser.add_argument('output', help='path of the JSON snapshot to write')
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as td:
        snap = build_snapshot(pathlib.Path(td))

    out = os.path.abspath(args.output)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, 'w') as fid:
        json.dump(snap, fid, indent=1, sort_keys=True)
        fid.write('\n')
    print('wrote {}'.format(args.output))


if __name__ == '__main__':
    main()
