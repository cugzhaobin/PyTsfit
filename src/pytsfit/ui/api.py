#!/usr/bin/env python
# ----------------------------------------------------------
# Thin wrappers around the existing pytsfit API for the UI.
#
# Everything here delegates to the core modules (data / models /
# tsfitting / output / PyTsfit); no fitting or plotting logic is
# re-implemented. The wrappers only marshal the UI's inputs into the
# shapes the core functions expect, and turn core outputs into
# DataFrames / plotly Figures for display.
# ----------------------------------------------------------
import glob, io, os, re
import numpy as np
import pandas as pd

from pytsfit.data import posData, neuData
from pytsfit.models import build_param_dict
from pytsfit.tsfitting import tsfitting
from pytsfit.output import plot_obs_mod
from pytsfit.output import (output_velo, output_eqoffset, output_break,
                            output_postseismic_disp, output_postseismic_ts)

import plotly.graph_objects as go
from plotly.subplots import make_subplots

COMPONENTS = ('N', 'E', 'U')
_COLORS = {'N': '#d62728', 'E': '#2ca02c', 'U': '#1f77b4'}


# --------------------------------------------------------------------------
# Data discovery / loading
# --------------------------------------------------------------------------
def list_sites(tsdir, tsformat, sitefile=''):
    '''
    Return the sorted site list for the UI.

    Uses the sitefile when it exists (authoritative, same as the CLI);
    otherwise falls back to scanning ``tsdir`` for ``*.{tsformat}`` files
    and taking the leading token of each filename (before the first ``.``
    or ``_``) as a candidate site.
    '''
    if sitefile and os.path.isfile(sitefile):
        sites = np.genfromtxt(sitefile, dtype=str)
        if sites.size == 1:
            sites = [str(sites)]
        return sorted({str(s) for s in sites})
    files = sorted(glob.glob(os.path.join(tsdir, '*.{}'.format(tsformat))))
    sites = sorted({re.split(r'[._]', os.path.basename(f))[0] for f in files})
    return sites


def load_site(tsdir, tsformat, site):
    '''
    Load one site's time series. Mirrors the CLI's file discovery:
    the first ``{tsdir}/{site}*.{tsformat}`` match wins.
    '''
    posfiles = sorted(glob.glob('{}/{}*.{}'.format(tsdir, site, tsformat)))
    if not posfiles:
        raise FileNotFoundError(
            'no {}.{} time series found for site {} under {}'
            .format(site, tsformat, site, tsdir))
    posfile = posfiles[0]
    if tsformat == 'pos':
        return posData(posfile)
    if tsformat == 'neu':
        return neuData(posfile)
    raise ValueError('unknown tsformat: {} (expected "pos" or "neu")'.format(tsformat))


def station_metadata(tsdir, tsformat, sites):
    '''
    Load lon/lat for every site so the station-distribution map can be
    drawn. Returns a DataFrame with columns site/lon/lat; sites that fail
    to load are skipped (their name is kept with NaN coordinates).
    '''
    rows = []
    for site in sites:
        try:
            data = load_site(tsdir, tsformat, site)
            rows.append({'site': site, 'lon': float(data.lon), 'lat': float(data.lat)})
        except Exception:
            rows.append({'site': site, 'lon': np.nan, 'lat': np.nan})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Fitting (reuses build_param_dict + tsfitting, like do_pytsfit)
# --------------------------------------------------------------------------
def run_fit(dict_input, dict_param, fit_opts, site):
    '''
    Fit all three components of ``site``. Mirrors the N/E/U block of
    do_pytsfit: each component gets its own eqlist/eqpostlist depending on
    the ``eqoffset_ne/eqoffset_up`` / ``eqpost_ne/eqpost_up`` switches.

    Returns a dict {'N': nrun, 'E': erun, 'U': urun} of fitted tsfitting
    instances, plus the loaded data as {'data': data}.
    '''
    tsdir   = dict_input['tsdir']
    tsformat = dict_input['tsformat']
    timespan = dict_input.get('timespan', [-np.inf, np.inf])
    data = load_site(tsdir, tsformat, site)

    base_param = build_param_dict(dict_param, dict_input['eqfile'],
                                  dict_input.get('prior_velfile', ''),
                                  dict_input.get('prior_offsetfile', ''),
                                  dict_input.get('prior_periodfile', ''))

    runs = {}
    specs = [('N', data.N, data.SN, 'eqoffset_ne', 'eqpost_ne'),
             ('E', data.E, data.SE, 'eqoffset_ne', 'eqpost_ne'),
             ('U', data.U, data.SU, 'eqoffset_up', 'eqpost_up')]
    for comp, obs, sig, eqkey, eqpostkey in specs:
        pd_ = dict(base_param)
        pd_['eqlist']     = base_param['eqlist'] if dict_param[eqkey] else []
        pd_['eqpostlist'] = base_param['eqpostlist'] if dict_param[eqpostkey] else []
        run = tsfitting(data.site, data.lon, data.lat, data.decyr, obs, sig,
                        pd_, comp, timespan, fit_opts=fit_opts)
        run.doFitting()
        runs[comp] = run
    return {'data': data, 'runs': runs}


# --------------------------------------------------------------------------
# Results as tables
# --------------------------------------------------------------------------
def param_table(runs):
    '''
    One row per estimated parameter: component, parameter name (from
    ``flag2``), value and standard error (sqrt of the covariance diagonal).
    '''
    rows = []
    for comp in COMPONENTS:
        run = runs.get(comp)
        if run is None or not hasattr(run, 'param') or len(run.param) == 0:
            continue
        if run.cov is not None and len(run.param) == len(run.flag2):
            stderr = np.sqrt(np.maximum(np.diag(run.cov), 0.0))
        else:
            stderr = np.full(len(run.param), np.nan)
        for name, value, err in zip(run.flag2, run.param, stderr):
            rows.append({'Component': comp, 'Parameter': name,
                         'Value': value, 'StdErr': err})
    return pd.DataFrame(rows)


def fit_summary(runs):
    '''
    Per-component quality summary: WRMS, NRMS, sigma scale, edited count.
    '''
    rows = []
    for comp in COMPONENTS:
        run = runs.get(comp)
        if run is None or not hasattr(run, 'wrms'):
            continue
        rows.append({'Component': comp,
                     'WRMS (mm)': run.wrms,
                     'NRMS': getattr(run, 'nrms', np.nan),
                     'sigma_scale': getattr(run, 'sig_scale', np.nan),
                     'n_used': int(np.sum(run.good)) if hasattr(run, 'good') else len(run.t),
                     'n_edited': getattr(run, 'nedit', 0)})
    return pd.DataFrame(rows)


def run_fit_batch(dict_input, dict_param, fit_opts, sites, progress_cb=None):
    '''
    Fit every site in ``sites`` by calling run_fit per site.

    ``progress_cb(i, n)`` is invoked before each site is fitted (i = 1-based
    index, n = total), so the UI can drive st.progress. Sites that raise are
    captured instead of aborting the batch.

    Returns {'results': {site: run_fit-output}, 'failures': [(site, error)]}.
    '''
    results, failures = {}, []
    n = len(sites)
    for i, site in enumerate(sites, start=1):
        if progress_cb is not None:
            progress_cb(i, n)
        try:
            results[site] = run_fit(dict_input, dict_param, fit_opts, site)
        except Exception as exc:  # keep fitting the remaining sites
            failures.append((site, str(exc)))
    return {'results': results, 'failures': failures}


def batch_velocity_table(results):
    '''
    One row per fitted site: lon/lat plus N/E/U secular velocities (mm/yr)
    with standard errors and WRMS, extracted from the VELOCITY parameter of
    each component (same source as output_velo). Sites whose component fit
    produced no velocity are kept with NaN for that component.
    '''
    rows = []
    for site, res in results.items():
        runs = res['runs']
        row = {'site': site, 'lon': runs['N'].lon, 'lat': runs['N'].lat}
        for comp in COMPONENTS:
            run = runs.get(comp)
            has_vel = (run is not None and hasattr(run, 'flag2')
                       and 'VELOCITY' in run.flag2)
            if has_vel:
                idx = np.where(run.flag2 == 'VELOCITY')[0][0]
                vel = float(run.param[idx])
                err = float(np.sqrt(np.maximum(np.diag(run.cov), 0.0))[idx]) \
                    if run.cov is not None else np.nan
                wrms = float(getattr(run, 'wrms', np.nan))
            else:
                vel = err = wrms = np.nan
            row['V{}'.format(comp)] = vel
            row['S{}'.format(comp)] = err
            row['W{}'.format(comp)] = wrms
        rows.append(row)
    return pd.DataFrame(rows).set_index('site')


def batch_output_velo(results, fmt='DETAIL'):
    '''
    Render the CLI's velocity table for every fitted site into one string by
    calling output_velo per site (same GMT/IOS3D/GLOBK/DETAIL format).
    '''
    buf = io.StringIO()
    for site in sorted(results):
        runs = results[site]['runs']
        try:
            output_velo(runs['N'], runs['E'], runs['U'], fid=buf, fmt=fmt)
        except Exception:
            continue
    return buf.getvalue()


def outputs_to_csv(runs, dict_param, dict_output, fmt_velo='DETAIL'):
    '''
    Render the same outputs the CLI writes to files (velocity, eq offsets,
    breaks, postseismic) into in-memory strings for the download buttons.
    Returns a dict name -> CSV text; entries whose model term is absent
    come back empty.
    '''
    nrun, erun, urun = runs['N'], runs['E'], runs['U']
    out = {}

    buf = io.StringIO()
    output_velo(nrun, erun, urun, fid=buf, fmt=fmt_velo)
    out['velfile'] = buf.getvalue()

    if dict_output.get('eqoffset'):
        buf = io.StringIO()
        output_eqoffset(nrun, erun, urun, fid=buf, fmt='GMT3D')
        out['eqoffset'] = buf.getvalue()
    if dict_output.get('break'):
        buf = io.StringIO()
        output_break(nrun, erun, urun, fid=buf)
        out['break'] = buf.getvalue()
    if dict_output.get('eqpostdisp'):
        buf = io.StringIO()
        if dict_output.get('eqpost_tspan'):
            output_postseismic_disp(nrun, erun, urun,
                                    list(dict_output['eqpost_tspan']), fid=buf)
        out['eqpostdisp'] = buf.getvalue()
    if dict_output.get('eqpostts'):
        buf = io.StringIO()
        if dict_output.get('eqpost_tspan'):
            output_postseismic_ts(nrun, erun, urun,
                                  list(dict_output['eqpost_tspan']))
        out['eqpostts'] = buf.getvalue()
    return out


# --------------------------------------------------------------------------
# Plotly figures (all rendering goes through plotly)
# --------------------------------------------------------------------------
def make_raw_ts_figure(data, component=None):
    '''
    Three-panel (N/E/U) plot of the observed time series with error bars.
    '''
    fig = make_subplots(rows=3, cols=1, shared_xaxes=True,
                        subplot_titles=('North (mm)', 'East (mm)', 'Vertical (mm)'))
    specs = [('N', data.N, data.SN), ('E', data.E, data.SE), ('U', data.U, data.SU)]
    for row, (comp, obs, sig) in enumerate(specs, start=1):
        fig.add_trace(go.Scatter(
            x=data.decyr, y=obs, mode='markers', name=comp,
            marker={'color': _COLORS[comp], 'size': 4},
            error_y={'type': 'data', 'array': sig, 'visible': True,
                     'thickness': 0.8, 'width': 0.5, 'color': 'rgba(0,0,0,0.35)'},
        ), row=row, col=1)
    fig.update_xaxes(title_text='Time (year)', row=3, col=1)
    fig.update_layout(height=1000, margin={'t': 60}, showlegend=False,
                      title='Raw time series')
    return fig


def make_obsmod_figure(runs, plot_dict):
    '''
    Observed + modeled N/E/U figure, reusing the existing plot_obs_mod
    plotly branch (outfile=False returns the Figure without writing disk).
    '''
    plot_dict = dict(plot_dict)
    plot_dict['figformat'] = 'html'
    plot_dict['showfig'] = False
    nrun, erun, urun = runs['N'], runs['E'], runs['U']
    return plot_obs_mod(nrun, erun, urun,
                        nrun.param, erun.param, urun.param,
                        plot_dict, outfile=False)


def make_residual_figure(runs):
    '''
    Post-fit residuals (obs - model) for each component, with a zero line.
    '''
    fig = make_subplots(rows=3, cols=1, shared_xaxes=True,
                        subplot_titles=('N residual (mm)', 'E residual (mm)',
                                        'U residual (mm)'))
    for row, comp in enumerate(COMPONENTS, start=1):
        run = runs[comp]
        if not hasattr(run, 'res') or len(run.res) == 0:
            continue
        fig.add_trace(go.Scatter(
            x=run.t, y=run.res, mode='markers', name=comp,
            marker={'color': _COLORS[comp], 'size': 4},
        ), row=row, col=1)
        fig.add_hline(y=0, line_dash='dot', line_color='black', row=row, col=1)
    fig.update_xaxes(title_text='Time (year)', row=3, col=1)
    fig.update_layout(height=800, margin={'t': 60}, showlegend=False,
                      title='Residuals')
    return fig


def make_sitemap(meta):
    '''
    Station-distribution map (go.Scattergeo, no mapbox token required).
    Sites with missing coordinates are dropped.
    '''
    meta = meta.dropna(subset=['lon', 'lat'])
    if meta.empty:
        return go.Figure()
    lon, lat = meta['lon'].values, meta['lat'].values
    pad = max(0.5, (max(lon) - min(lon)) * 0.1, (max(lat) - min(lat)) * 0.1)

    fig = go.Figure(go.Scattergeo(
        lon=lon, lat=lat, mode='markers+text',
        text=meta['site'].values, textposition='top center',
        textfont={'size': 10},
        marker={'size': 8, 'color': '#d62728', 'line': {'width': 1, 'color': 'black'}},
        customdata=meta['site'].values,
        hovertemplate='%{customdata}<br>lon %{lon:.3f} lat %{lat:.3f}<extra></extra>',
    ))
    fig.update_layout(
        height=650, margin={'t': 60, 'b': 0, 'l': 0, 'r': 0},
        title='Station distribution',
        geo={'projection_type': 'equirectangular',
             'showland': True, 'landcolor': 'rgb(235,235,235)',
             'showcoastlines': True, 'coastlinecolor': 'rgb(160,160,160)',
             'showcountries': True, 'countrycolor': 'rgb(190,190,190)',
             'lonaxis': {'range': [min(lon) - pad, max(lon) + pad]},
             'lataxis': {'range': [min(lat) - pad, max(lat) + pad]}},
    )
    return fig
