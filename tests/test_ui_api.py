'''
Tests for the UI layer (pytsfit.ui.api): the thin wrappers that marshal
Streamlit inputs into the existing pytsfit API and back. Everything here
delegates to the core modules, so these tests guard the wiring, not the
fitting/plotting algorithms themselves (which the characterization tests
cover).
'''
import os

import numpy as np
import pytest

from fixtures import write_pos, write_eq_rename

from pytsfit.qualitycontrol import DEFAULT_FIT_OPTS
from pytsfit.ui import api

import plotly.graph_objects as go


@pytest.fixture()
def tsdir(tmp_path):
    d = tmp_path / 'pos'
    d.mkdir()
    write_pos(str(d / 'TEST.pos'), site='TEST')
    write_pos(str(d / 'TESTB.pos'), site='TESTB')
    return str(d)


@pytest.fixture()
def dict_input(tsdir, tmp_path):
    return {
        'tsdir': tsdir, 'tsformat': 'pos', 'sitefile': '',
        'eqfile': str(write_eq_rename(str(tmp_path / 'eq_rename'))),
        'prior_velfile': '', 'prior_offsetfile': '', 'prior_periodfile': '',
        'timespan': [1998, 2024],
    }


def _full_param():
    return {'constant': True, 'linear': True, 'annual': True, 'semiannual': True,
            'break': True, 'eqoffset_ne': True, 'eqoffset_up': True,
            'eqpost_ne': True, 'eqpost_up': True}


def test_list_sites_from_directory(tsdir):
    sites = api.list_sites(tsdir, 'pos', '')
    assert set(sites) == {'TEST', 'TESTB'}


def test_list_sites_from_sitefile(tsdir, tmp_path):
    sitefile = tmp_path / 'sites.lst'
    sitefile.write_text('TEST\n')
    assert api.list_sites(tsdir, 'pos', str(sitefile)) == ['TEST']


def test_load_site(tsdir):
    data = api.load_site(tsdir, 'pos', 'TEST')
    assert data.site == 'TEST'
    assert len(data.decyr) > 100
    assert data.lon == pytest.approx(103.2, abs=1e-3)


def test_station_metadata(tsdir):
    meta = api.station_metadata(tsdir, 'pos', ['TEST', 'TESTB'])
    assert set(meta['site']) == {'TEST', 'TESTB'}
    assert meta['lon'].notna().all()


def test_run_fit_recovers_known_parameters(dict_input):
    res = api.run_fit(dict_input, _full_param(), dict(DEFAULT_FIT_OPTS), 'TEST')
    run = res['runs']['N']
    assert 'VELOCITY' in run.flag
    # Synthetic series: N velocity = -3.5 mm/yr, coseismic step = +18 mm.
    assert run.param[1] == pytest.approx(-3.5, abs=0.2)
    assert run.param[2] == pytest.approx(18.0, abs=1.0)
    assert run.wrms < 2.0
    assert set(res['runs']) == {'N', 'E', 'U'}


def test_param_table_and_summary(dict_input):
    res = api.run_fit(dict_input, _full_param(), dict(DEFAULT_FIT_OPTS), 'TEST')
    table = api.param_table(res['runs'])
    assert set(table['Component']) == {'N', 'E', 'U'}
    assert {'Component', 'Parameter', 'Value', 'StdErr'} <= set(table.columns)
    assert not table['StdErr'].isna().all()

    summary = api.fit_summary(res['runs'])
    assert len(summary) == 3
    assert {'WRMS (mm)', 'NRMS', 'n_used'} <= set(summary.columns)


def test_figures_returned_not_written(dict_input):
    res = api.run_fit(dict_input, _full_param(), dict(DEFAULT_FIT_OPTS), 'TEST')
    runs = res['runs']
    data = res['data']
    plot_dict = {'detrend': True, 'debreak': True, 'deeqoffset': False,
                 'depost': False, 'deseason': False}

    fig = api.make_obsmod_figure(runs, plot_dict)
    assert isinstance(fig, go.Figure)
    # outfile=False path must not write the legacy {site}.html to disk.
    assert not os.path.exists(os.path.join(dict_input['tsdir'], 'TEST.html'))

    assert isinstance(api.make_residual_figure(runs), go.Figure)
    assert isinstance(api.make_raw_ts_figure(data), go.Figure)
    meta = api.station_metadata(dict_input['tsdir'], 'pos', ['TEST'])
    assert isinstance(api.make_sitemap(meta), go.Figure)


def test_run_fit_batch_fits_all_sites(dict_input):
    progress = []
    batch = api.run_fit_batch(dict_input, _full_param(), dict(DEFAULT_FIT_OPTS),
                              ['TEST', 'TESTB'],
                              progress_cb=lambda i, n: progress.append((i, n)))
    assert set(batch['results']) == {'TEST', 'TESTB'}
    assert batch['failures'] == []
    assert progress == [(1, 2), (2, 2)]
    # N velocity of the synthetic series is -3.5 mm/yr.
    assert batch['results']['TEST']['runs']['N'].param[1] == pytest.approx(-3.5, abs=0.2)


def test_run_fit_batch_keeps_going_on_failure(dict_input):
    batch = api.run_fit_batch(dict_input, _full_param(), dict(DEFAULT_FIT_OPTS),
                              ['TEST', 'NOPE'])
    assert set(batch['results']) == {'TEST'}
    assert [s for s, _ in batch['failures']] == ['NOPE']


def test_batch_velocity_table(dict_input):
    batch = api.run_fit_batch(dict_input, _full_param(), dict(DEFAULT_FIT_OPTS),
                              ['TEST', 'TESTB'])
    table = api.batch_velocity_table(batch['results'])
    assert list(table.index) == ['TEST', 'TESTB']
    assert {'VN', 'SN', 'VE', 'SE', 'VU', 'SU', 'WN', 'WE', 'WU'} <= set(table.columns)
    # N velocity = -3.5 mm/yr with a small standard error.
    assert table.loc['TEST', 'VN'] == pytest.approx(-3.5, abs=0.2)
    assert table.loc['TEST', 'SN'] < 1.0


def test_batch_output_velo(dict_input):
    batch = api.run_fit_batch(dict_input, _full_param(), dict(DEFAULT_FIT_OPTS),
                              ['TEST', 'TESTB'])
    text = api.batch_output_velo(batch['results'], fmt='DETAIL')
    assert 'TEST' in text and 'TESTB' in text
    assert len(text.strip().splitlines()) == 2


def test_outputs_to_csv(dict_input):
    dp = _full_param()
    res = api.run_fit(dict_input, dp, dict(DEFAULT_FIT_OPTS), 'TEST')
    dout = {'velfile': True, 'eqoffset': True, 'break': True,
            'eqpostdisp': False, 'eqpostts': False, 'eqpost_tspan': [2015, 2023]}
    csvs = api.outputs_to_csv(res['runs'], dp, dout)
    assert 'TEST' in csvs['velfile']
    assert csvs['eqoffset'].strip() != ''
    # break lines carry E/N offsets + epoch (no site name in this format)
    assert '2012.208' in csvs['break'] and '2016.669' in csvs['break']
