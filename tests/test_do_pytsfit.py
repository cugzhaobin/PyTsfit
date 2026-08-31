'''
Tests for the YAML / command-line override logic in do_pytsfit.py.

The script reads its main parameters from a YAML config file. Command-line
options override those parameters only when the user actually supplies them,
so the YAML file remains the single source of defaults. These tests pin that
behaviour by exercising build_parser() and apply_cli_overrides() directly.
'''
import argparse

import pytest

from pytsfit.scripts import do_pytsfit


def _base_cfg():
    '''A config dict shaped like config.yaml.'''
    return {
        'dict_input': {'eqfile': '', 'prior_velfile': '', 'prior_offsetfile': '',
                       'prior_periodfile': '', 'sitefile': '', 'timespan': [1998, 2024],
                       'tsdir': '../pos/', 'tsformat': 'pos'},
        'dict_param': {'constant': True, 'linear': True, 'annual': True,
                       'semiannual': False, 'break': False, 'eqoffset_ne': False,
                       'eqoffset_up': False, 'eqpost_ne': False, 'eqpost_up': False},
        'dict_plot': {'detrend': True, 'debreak': True, 'deeqoffset': False,
                      'depost': False, 'deseason': False, 'showfig': False,
                      'figformat': 'jpg'},
        'dict_output': {'tsfig': False, 'param': False, 'obsmod': False,
                        'eqpostts': False, 'velfile': '', 'eqoffset': '',
                        'break': '', 'eqpostdisp': '', 'eqpost_tspan': [2015, 2023]},
        'dict_fit': {'sigma_scale': 'nrms', 'min_rsig': 30, 'min_sigscale': 10,
                     'outlier': False, 'nsigma': 4.0, 'outlier_scale': 'mad',
                     'max_iter': 10, 'restore_factor': 0.9, 'max_sigma': 0.0},
    }


def _parse(argv):
    return do_pytsfit.build_parser().parse_args(argv)


# ---------------------------------------------------------------------------
# build_parser: every YAML-mirroring option defaults to None, so the absence of
# a flag is distinguishable from an explicit value.
# ---------------------------------------------------------------------------
def test_parser_defaults_are_none():
    args = _parse(['--cfgfile', 'x.yaml'])
    for key, _ in (do_pytsfit._BOOL_KEYS + do_pytsfit._STR_KEYS
                   + do_pytsfit._LIST_KEYS + do_pytsfit._FLOAT_KEYS
                   + do_pytsfit._INT_KEYS):
        assert getattr(args, key) is None, key
    assert args.sitelist is None


def test_parser_lowercases_bool_values():
    args = _parse(['--cfgfile', 'x.yaml', '--annual', 'TRUE', '--showfig', 'False'])
    assert args.annual == 'true'
    assert args.showfig == 'false'


def test_parser_rejects_invalid_bool(tmp_path):
    with pytest.raises(SystemExit):
        _parse(['--cfgfile', 'x.yaml', '--annual', 'yes'])


def test_parser_parses_list_flags_as_float_pairs():
    args = _parse(['--cfgfile', 'x.yaml', '--timespan', '2000', '2021',
                   '--eqpost-tspan', '2015.5', '2022.5'])
    assert args.timespan == [2000.0, 2021.0]
    assert args.eqpost_tspan == [2015.5, 2022.5]


# ---------------------------------------------------------------------------
# apply_cli_overrides: non-None options overwrite the YAML value, None leaves it.
# ---------------------------------------------------------------------------
def test_no_cli_options_leaves_config_untouched():
    cfg = _base_cfg()
    args = _parse(['--cfgfile', 'x.yaml'])
    do_pytsfit.apply_cli_overrides(args, cfg)
    assert cfg['dict_param']['annual'] is True
    assert cfg['dict_param']['semiannual'] is False
    assert cfg['dict_plot']['showfig'] is False
    assert cfg['dict_input']['timespan'] == [1998, 2024]
    assert cfg['dict_input']['tsdir'] == '../pos/'
    assert cfg['dict_input']['prior_velfile'] == ''
    assert cfg['dict_output']['tsfig'] is False
    assert cfg['dict_output']['velfile'] == ''
    assert cfg['dict_output']['eqpost_tspan'] == [2015, 2023]


def test_cli_overrides_bool_values():
    cfg = _base_cfg()
    args = _parse(['--cfgfile', 'x.yaml', '--annual', 'false',
                   '--semiannual', 'true', '--showfig', 'true',
                   '--tsfig', 'true', '--param', 'true'])
    do_pytsfit.apply_cli_overrides(args, cfg)
    assert cfg['dict_param']['annual'] is False
    assert cfg['dict_param']['semiannual'] is True
    assert cfg['dict_plot']['showfig'] is True
    assert cfg['dict_output']['tsfig'] is True
    assert cfg['dict_output']['param'] is True


def test_cli_overrides_input_strings_and_timespan():
    cfg = _base_cfg()
    args = _parse(['--cfgfile', 'x.yaml', '--tsdir', '../data/pos/',
                   '--tsformat', 'neu', '--timespan', '2005', '2019'])
    do_pytsfit.apply_cli_overrides(args, cfg)
    assert cfg['dict_input']['tsdir'] == '../data/pos/'
    assert cfg['dict_input']['tsformat'] == 'neu'
    assert cfg['dict_input']['timespan'] == [2005.0, 2019.0]


def test_cli_overrides_prior_files():
    cfg = _base_cfg()
    args = _parse(['--cfgfile', 'x.yaml', '--prior-velfile', 'velo.vel',
                   '--prior-offsetfile', 'offset.dat',
                   '--prior-periodfile', 'season.dat'])
    do_pytsfit.apply_cli_overrides(args, cfg)
    assert cfg['dict_input']['prior_velfile'] == 'velo.vel'
    assert cfg['dict_input']['prior_offsetfile'] == 'offset.dat'
    assert cfg['dict_input']['prior_periodfile'] == 'season.dat'
    # the output velocity file is a distinct key and stays untouched
    assert cfg['dict_output']['velfile'] == ''


def test_cli_overrides_output_paths_and_tspan():
    cfg = _base_cfg()
    args = _parse(['--cfgfile', 'x.yaml', '--velfile', 'velo.gmtvec',
                   '--break', 'break.out', '--eqpost-tspan', '2016', '2021'])
    do_pytsfit.apply_cli_overrides(args, cfg)
    assert cfg['dict_output']['velfile'] == 'velo.gmtvec'
    assert cfg['dict_output']['break'] == 'break.out'
    assert cfg['dict_output']['eqpost_tspan'] == [2016.0, 2021.0]


def test_partial_override_touches_only_given_options():
    cfg = _base_cfg()
    args = _parse(['--cfgfile', 'x.yaml', '--showfig', 'true'])
    do_pytsfit.apply_cli_overrides(args, cfg)
    assert cfg['dict_plot']['showfig'] is True
    # everything else keeps its YAML value
    assert cfg['dict_param']['annual'] is True
    assert cfg['dict_param']['semiannual'] is False
    assert cfg['dict_input']['timespan'] == [1998, 2024]
    assert cfg['dict_output']['velfile'] == ''


# ---------------------------------------------------------------------------
# dict_fit: the new fitting/quality-control section (added for the
# realistic-sigma and outlier features).
# ---------------------------------------------------------------------------
def test_parser_parses_fit_numeric_flags():
    args = _parse(['--cfgfile', 'x.yaml', '--nsigma', '5.5',
                   '--restore-factor', '0.8', '--max-sigma', '20.0',
                   '--max-iter', '7', '--min-rsig', '40', '--min-sigscale', '12'])
    assert args.nsigma == 5.5
    assert args.restore_factor == 0.8
    assert args.max_sigma == 20.0
    assert args.max_iter == 7
    assert args.min_rsig == 40
    assert args.min_sigscale == 12


def test_cli_overrides_fit_section():
    cfg = _base_cfg()
    args = _parse(['--cfgfile', 'x.yaml', '--sigma-scale', 'realistic',
                   '--outlier-scale', 'sigma', '--outlier', 'true',
                   '--nsigma', '3.0', '--max-iter', '5'])
    do_pytsfit.apply_cli_overrides(args, cfg)
    assert cfg['dict_fit']['sigma_scale'] == 'realistic'
    assert cfg['dict_fit']['outlier_scale'] == 'sigma'
    assert cfg['dict_fit']['outlier'] is True
    assert cfg['dict_fit']['nsigma'] == 3.0
    assert cfg['dict_fit']['max_iter'] == 5


def test_no_cli_fit_options_leaves_fit_section_untouched():
    cfg = _base_cfg()
    args = _parse(['--cfgfile', 'x.yaml'])
    do_pytsfit.apply_cli_overrides(args, cfg)
    assert cfg['dict_fit']['sigma_scale'] == 'nrms'
    assert cfg['dict_fit']['outlier'] is False
    assert cfg['dict_fit']['nsigma'] == 4.0
