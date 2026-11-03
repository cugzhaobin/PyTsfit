#!/usr/bin/env python
import glob, yaml, argparse, os, logging
import numpy as np
from pytsfit.data import posData, neuData
from pytsfit.models import build_param_dict
from pytsfit.tsfitting import tsfitting
from pytsfit.output import (plot_obs_mod, output_param, output_obs_mod,
                            output_velo, output_eqoffset, output_break,
                            output_postseismic_disp, output_postseismic_ts,
                            pair_postseismic_outputs, open_postdisp_file,
                            resolve_displacement_spans,
                            resolve_postseismic_ts_span)

# Header lines describing the columns of each output file. These are written
# as '#' comment lines so they are ignored by GMT and np.genfromtxt. The
# postseismic displacement header is built per window by output.postdisp_header().
HEADER_VELO      = '# Lon(deg)  Lat(deg)  Ve(mm/yr)  Vn(mm/yr)  Ve(mm/yr)  Vn(mm/yr)  Se(mm/yr)  Sn(mm/yr)  Corr  Vu(mm/yr)  Vu(mm/yr)  Su(mm/yr)  Site  #  tstart(yr)  tend(yr)  tspan(yr)  nobs\n'
HEADER_EQOFFSET  = '# Lon(deg)  Lat(deg)  E_off(mm)  N_off(mm)  U_off(mm)  Se(mm)  Sn(mm)  Su(mm)  Site  #  EQ_code\n'

def write_header_if_empty(fid, fname, header):
    '''
    Write a comment header line to the output file if it is empty.
    The output files are opened in append mode, so this avoids writing a
    duplicate header on subsequent runs of the same script.
    '''
    if os.path.getsize(fname) == 0:
        fid.write(header)

# YAML keys exposed as command-line overrides, grouped by value type. Each
# entry is ``(key, section)``; the CLI flag is ``--<key>`` with underscores
# turned into dashes (e.g. ``eqpostdsp_tspan`` -> ``--eqpostdsp-tspan``), and its
# ``dest`` is ``key`` so it can be read back with ``getattr(args, key)``.
_BOOL_KEYS = [
    ('annual',     'dict_param'),
    ('semiannual', 'dict_param'),
    ('showfig',    'dict_plot'),
    ('tsfig',      'dict_output'),
    ('param',      'dict_output'),
    ('obsmod',     'dict_output'),
    ('eqpostts',   'dict_output'),
    ('outlier',    'dict_fit'),
]

_STR_KEYS = [
    ('tsdir',           'dict_input'),
    ('tsformat',        'dict_input'),
    ('prior_velfile',   'dict_input'),
    ('prior_offsetfile','dict_input'),
    ('prior_periodfile','dict_input'),
    ('velfile',         'dict_output'),
    ('eqoffset',        'dict_output'),
    ('break',           'dict_output'),
    ('sigma_scale',     'dict_fit'),
    ('outlier_scale',   'dict_fit'),
]

# Keys that accept one or more values. `eqpostdisp` pairs with
# `eqpostdsp_tspan` positionally, so both support multiple values on the
# command line.
_MULTI_STR_KEYS = [
    ('eqpostdisp', 'dict_output'),
]

_LIST_KEYS = [
    ('timespan',      'dict_input'),
    ('eqpostts_tspan', 'dict_output'),
]

_MULTI_LIST_KEYS = [
    ('eqpostdsp_tspan', 'dict_output'),
]

# Additional command-line spellings kept for backward compatibility. The
# canonical option is built from the key; these aliases point at the same
# destination so an older command line keeps working.
_OPTION_ALIASES = {
    'eqpostdsp_tspan': ['--eqpost-tspan'],
}

# Numeric fit options exposed as command-line overrides. Kept separate from the
# string keys because their values must be parsed as floats/ints, not left as
# strings.
_FLOAT_KEYS = [
    ('nsigma',         'dict_fit'),
    ('restore_factor', 'dict_fit'),
    ('max_sigma',      'dict_fit'),
]

_INT_KEYS = [
    ('max_iter',    'dict_fit'),
    ('min_rsig',    'dict_fit'),
    ('min_sigscale','dict_fit'),
]

def build_parser():
    '''
    Build the command-line parser.

    Every option that mirrors a YAML parameter defaults to ``None``. A ``None``
    means "not given on the command line", so the value read from the YAML file
    is used unchanged; a non-``None`` value overrides that YAML parameter. This
    keeps the YAML file as the single source of defaults.
    '''
    parser = argparse.ArgumentParser(description="Position time series fitting.")
    parser.add_argument('--cfgfile', type=str, required=True, help='configure file in YAML format')
    parser.add_argument('--sitelist', type=str, required=False, help='This will overwrite the sitefile in the configure file.', nargs='+')
    for key, section in _BOOL_KEYS:
        parser.add_argument('--' + key.replace('_', '-'),
                            type=str.lower, choices=['true', 'false'],
                            help='Override {}.{} ("true" or "false").'.format(section, key))
    for key, section in _STR_KEYS:
        parser.add_argument('--' + key.replace('_', '-'), type=str,
                            help='Override {}.{}.'.format(section, key))
    for key, section in _MULTI_STR_KEYS:
        parser.add_argument('--' + key.replace('_', '-'), type=str, nargs='+',
                            help='Override {}.{} (one or more values).'
                                 .format(section, key))
    for key, section in _LIST_KEYS:
        parser.add_argument('--' + key.replace('_', '-'), type=float, nargs=2,
                            help='Override {}.{} (two values).'.format(section, key))
    for key, section in _MULTI_LIST_KEYS:
        option = '--' + key.replace('_', '-')
        aliases = [option] + _OPTION_ALIASES.get(key, [])
        parser.add_argument(*aliases, dest=key, type=float, nargs='+',
                            help='Override {}.{} (t0 t1 [t0 t1 ...]).'
                                 .format(section, key))
    for key, section in _FLOAT_KEYS:
        parser.add_argument('--' + key.replace('_', '-'), type=float,
                            help='Override {}.{}.'.format(section, key))
    for key, section in _INT_KEYS:
        parser.add_argument('--' + key.replace('_', '-'), type=int,
                            help='Override {}.{}.'.format(section, key))
    return parser

def apply_cli_overrides(args, cfg):
    '''
    Overlay command-line options onto the YAML config.

    Only options the user actually supplied (non-``None``) overwrite the
    corresponding YAML value; everything else is left untouched. The config is
    mutated in place and returned for convenience.
    '''
    for key, section in _BOOL_KEYS:
        value = getattr(args, key)
        if value is not None:
            cfg[section][key] = (value == 'true')
    for key, section in _STR_KEYS:
        value = getattr(args, key)
        if value is not None:
            cfg[section][key] = value
    for key, section in _MULTI_STR_KEYS:
        value = getattr(args, key)
        if value is not None:
            cfg[section][key] = list(value)
    for key, section in _LIST_KEYS:
        value = getattr(args, key)
        if value is not None:
            cfg[section][key] = list(value)
    for key, section in _MULTI_LIST_KEYS:
        value = getattr(args, key)
        if value is not None:
            if len(value) % 2 != 0:
                raise SystemExit(
                    '--{} expects an even number of values: '
                    't0 t1 [t0 t1 ...]'.format(key.replace('_', '-')))
            if len(value) == 2:
                # Keep the historical flat [t0, t1] form for one window.
                cfg[section][key] = [float(v) for v in value]
            else:
                cfg[section][key] = [[float(v) for v in value[i:i + 2]]
                                     for i in range(0, len(value), 2)]
    for key, section in _FLOAT_KEYS:
        value = getattr(args, key)
        if value is not None:
            cfg[section][key] = value
    for key, section in _INT_KEYS:
        value = getattr(args, key)
        if value is not None:
            cfg[section][key] = value
    return cfg

def main():
    logging.basicConfig(level=logging.INFO,
        format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
        datefmt="%d-%m-%Y %H:%M:%S")
    args = build_parser().parse_args()

    with open(args.cfgfile, 'r') as fid:
        lines = fid.read()
        cfg   = yaml.load(lines, Loader=yaml.FullLoader)

    # Command-line options override the YAML parameters, but only where the
    # user supplied them; otherwise the YAML value is kept.
    apply_cli_overrides(args, cfg)

    dict_input = cfg['dict_input']
    eqfile     = dict_input['eqfile']
    prior_velfile    = dict_input['prior_velfile']
    prior_offsetfile = dict_input['prior_offsetfile']
    prior_periodfile = dict_input['prior_periodfile']
    sitefile   = dict_input['sitefile']
    tsdir      = dict_input['tsdir']
    tsformat   = dict_input['tsformat']
    timespan   = dict_input['timespan']
    if os.path.isfile(sitefile) == True:
        sitelist   = np.genfromtxt(sitefile, dtype=str)
    else:
        sitelist   = np.array([])
    if sitelist.size == 1:
        sitelist = [str(sitelist)]
    if args.sitelist is not None:
        sitelist = args.sitelist

    dict_param = cfg['dict_param']
    dict_plot  = cfg['dict_plot']
    dict_output= cfg['dict_output']
    # Older config files may not have a dict_fit section; fall back to defaults.
    fit_opts   = cfg.get('dict_fit', {})
    param_dict = build_param_dict(dict_param, eqfile,
                                  prior_velfile, prior_offsetfile, prior_periodfile)

    fid_velo      = None
    fid_eq        = None
    fid_brk       = None
    if len(dict_output['velfile'])>0:
        fid_velo = open(dict_output['velfile'], 'a')
        write_header_if_empty(fid_velo, dict_output['velfile'], HEADER_VELO)
    if len(dict_output['eqoffset'])>0:
        fid_eq   = open(dict_output['eqoffset'], 'a')
        write_header_if_empty(fid_eq, dict_output['eqoffset'], HEADER_EQOFFSET)
    if len(dict_output['break'])>0:
        fid_brk  = open(dict_output['break'], 'a')

    # Pair each postseismic output file with its own integration window, so
    # several cumulative-displacement periods can be written side by side.
    # `eqpostdsp_tspan` is the current key; `eqpost_tspan` is accepted as an
    # older alias for the same setting.
    dsp_spans = resolve_displacement_spans(
        dict_output.get('eqpostdsp_tspan', []),
        dict_output.get('eqpost_tspan', []))
    postdisp_outputs = pair_postseismic_outputs(
        dict_output.get('eqpostdisp', ''), dsp_spans)
    fids_post_disp = [open_postdisp_file(fname, time_span)
                      for fname, time_span in postdisp_outputs]

    # The postseismic time-series output has its own window, independent of
    # the cumulative-displacement integration windows. Resolve it only when
    # that output is enabled, so a disabled output cannot warn or fail.
    postts_span = []
    if dict_output.get('eqpostts', False):
        postts_span = resolve_postseismic_ts_span(
            dict_output.get('eqpostts_tspan', []), dsp_spans)
        if not postts_span:
            logging.warning('eqpostts is enabled but no postseismic '
                            'time-series window is configured; skipping that '
                            'output.')

    # Optional 2-char earthquake code (eq_rename.eq) selecting the event whose
    # postseismic displacement is reported. Empty means "sum all EQDECAY terms".
    eqpost_eqcode = dict_output.get('eqpost_eqcode', '') or None

    for i in range(len(sitelist)):
        print(60*'-')
        posfiles = glob.glob('{}/{}*.{}'.format(tsdir, sitelist[i], tsformat))
        for j in range(len(posfiles)):
            logging.info('fitting time series for {}'.format(posfiles[j]))
            posfile = posfiles[j]
            if tsformat == 'pos':
                data    = posData(posfile)
            if tsformat == 'neu':
                data    = neuData(posfile)

            # North
            if dict_param['eqoffset_ne'] == False: param_dict['eqlist']=[]
            if dict_param['eqpost_ne'] == False: param_dict['eqpostlist']=[]
            nrun   = tsfitting(data.site, data.lon, data.lat, data.decyr, data.N, data.SN, param_dict, 'N', timespan, fit_opts=fit_opts)
            nparam = nrun.doFitting()
            print(nparam)

            # East
            if dict_param['eqoffset_ne'] == False: param_dict['eqlist']=[]
            if dict_param['eqpost_ne'] == False: param_dict['eqpostlist']=[]
            erun  = tsfitting(data.site, data.lon, data.lat, data.decyr, data.E, data.SE, param_dict, 'E', timespan, fit_opts=fit_opts)
            eparam = erun.doFitting()

            # Vertical
            if dict_param['eqoffset_up'] == False: param_dict['eqlist']=[]
            if dict_param['eqpost_up'] == False: param_dict['eqpostlist']=[]
            urun  = tsfitting(data.site, data.lon, data.lat, data.decyr, data.U, data.SU, param_dict, 'U', timespan, fit_opts=fit_opts)
            uparam = urun.doFitting()

            if len(nparam) == 0 or len(eparam) == 0 or len(uparam) == 0:
                continue

            if dict_output['tsfig'] == True:
                plot_obs_mod(nrun, erun, urun, nparam, eparam, uparam, dict_plot, nsigma=3, nwrms=3)
            if dict_output['param'] == True:
                output_param(nrun, erun, urun, nparam, eparam, uparam)
            if dict_output['obsmod'] == True:
                output_obs_mod(nrun, erun, urun, nparam, eparam, uparam, dict_plot)
            if len(dict_output['velfile']) > 0:
                output_velo(nrun, erun, urun, fid=fid_velo, fmt='DETAIL')
            if len(dict_output['eqoffset']) > 0:
                output_eqoffset(nrun, erun, urun, fid=fid_eq, fmt='GMT3D')
            if len(dict_output['break']) > 0:
                output_break(nrun, erun, urun, fid=fid_brk)
            for fid_post_disp, (_, post_tspan) in zip(fids_post_disp, postdisp_outputs):
                output_postseismic_disp(nrun, erun, urun, post_tspan,
                                        fid=fid_post_disp, eqcode=eqpost_eqcode)
            if dict_output['eqpostts'] == True and len(postts_span) != 0:
                output_postseismic_ts(nrun, erun, urun, postts_span,
                                      eqcode=eqpost_eqcode)
    if fid_velo      != None: fid_velo.close()
    if fid_eq        != None: fid_eq.close()
    if fid_brk       != None: fid_brk.close()
    for fid_post_disp in fids_post_disp:
        fid_post_disp.close()

if __name__ == '__main__':
    main()
