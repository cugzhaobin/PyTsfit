#!/usr/bin/env python
import glob, yaml, argparse, os, logging
import numpy as np
from pytsfit.PyTsfit import *

# Header lines describing the columns of each output file. These are written
# as '#' comment lines so they are ignored by GMT and np.genfromtxt.
HEADER_VELO      = '# Lon(deg)  Lat(deg)  Ve(mm/yr)  Vn(mm/yr)  Ve(mm/yr)  Vn(mm/yr)  Se(mm/yr)  Sn(mm/yr)  Corr  Vu(mm/yr)  Vu(mm/yr)  Su(mm/yr)  Site  #  tstart(yr)  tend(yr)  tspan(yr)  nobs\n'
HEADER_EQOFFSET  = '# Lon(deg)  Lat(deg)  E_off(mm)  N_off(mm)  U_off(mm)  Se(mm)  Sn(mm)  Su(mm)  Site  #  EQ_code\n'
HEADER_POSTDISP  = '# Lon(deg)  Lat(deg)  E_disp(mm)  N_disp(mm)  E_wrms(mm)  N_wrms(mm)  Corr  Site  U_disp(mm)  U_wrms(mm)\n'

def write_header_if_empty(fid, fname, header):
    '''
    Write a comment header line to the output file if it is empty.
    The output files are opened in append mode, so this avoids writing a
    duplicate header on subsequent runs of the same script.
    '''
    if os.path.getsize(fname) == 0:
        fid.write(header)

def main():
    logging.basicConfig(level=logging.INFO,
        format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
        datefmt="%d-%m-%Y %H:%M:%S")
    parser = argparse.ArgumentParser(description="Position time series fitting.")
    parser.add_argument('--cfgfile', type=str, required=True, help='configure file in YAML format')
    parser.add_argument('--sitelist', type=str, required=False, help='This will overwrite the sitefile in the configure file.', nargs='+')
    parser.add_argument('--showfig', type=str, required=False, help='Show figure.')
    parser.add_argument('--annual', type=str, required=False, help='Estimate annual term.')
    parser.add_argument('--semiannual', type=str, required=False, help='Estimate semiannual term.')
    args = parser.parse_args()

    with open(args.cfgfile, 'r') as fid:
        lines = fid.read()
        cfg   = yaml.load(lines, Loader=yaml.FullLoader)
    dict_input = cfg['dict_input']
    eqfile     = dict_input['eqfile']
    velfile    = dict_input['velfile']
    offsetfile = dict_input['offsetfile']
    periodfile = dict_input['periodfile']
    sitefile   = dict_input['sitefile']
    tsdir      = dict_input['tsdir']
    tsformat   = dict_input['tsformat']
    constraint = correction(velfile, offsetfile, periodfile)
    timespan   = dict_input['timespan']
    if os.path.isfile(sitefile) == True:
        sitelist   = np.genfromtxt(sitefile, dtype=str)
    else:
        sitelist   = np.array([])
    if sitelist.size == 1:
        sitelist = [str(sitelist)]
    if args.sitelist != None:
        sitelist = args.sitelist

    dict_param = cfg['dict_param']
    dict_plot  = cfg['dict_plot']
    dict_output= cfg['dict_output']
    param_dict = {
              'constant'  : dict_param['constant'],
              'linear'    : dict_param['linear'],
              'ANN'       : dict_param['annual'],
              'SANN'      : dict_param['semiannual'],
              'eqlist'    : [],
              'eqpostlist': [],
              'brklist'   : [],
              'correct'   : constraint}
    if args.annual is not None:
        param_dict['ANN'] = args.annual.lower() == 'true'
    if args.semiannual is not None:
        param_dict['SANN'] = args.semiannual.lower() == 'true'
    if dict_param['eqoffset_ne'] == True or dict_param['eqoffset_up'] == True:
        eq = eqcatalog(eqfile)
        param_dict['eqlist'] = eq.eqlist
    if dict_param['break'] == True:
        bk = breakcatalog(eqfile)
        param_dict['brklist'] = bk.breaklist
    if dict_param['eqpost_ne'] == True or dict_param['eqpost_up'] == True:
#       eq  = eqcatalog(eqfile)
        eqp = eqPostList(eqfile, eq)
        param_dict['eqpostlist'] = eqp.eqpostlist

    fid_velo      = None
    fid_eq        = None
    fid_brk       = None
    fid_post_disp = None
    if len(dict_output['velfile'])>0:
        fid_velo = open(dict_output['velfile'], 'a')
        write_header_if_empty(fid_velo, dict_output['velfile'], HEADER_VELO)
    if len(dict_output['eqoffset'])>0:
        fid_eq   = open(dict_output['eqoffset'], 'a')
        write_header_if_empty(fid_eq, dict_output['eqoffset'], HEADER_EQOFFSET)
    if len(dict_output['break'])>0:
        fid_brk  = open(dict_output['break'], 'a')
    if len(dict_output['eqpostdisp']) > 0:
        fid_post_disp = open(dict_output['eqpostdisp'], 'a')
        write_header_if_empty(fid_post_disp, dict_output['eqpostdisp'], HEADER_POSTDISP)

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
            nrun   = tsfitting(data.site, data.lon, data.lat, data.decyr, data.N, data.SN, param_dict, 'N', timespan)
            nparam = nrun.doFitting()
            print(nparam)
    	
            # East
            if dict_param['eqoffset_ne'] == False: param_dict['eqlist']=[]
            if dict_param['eqpost_ne'] == False: param_dict['eqpostlist']=[]
            erun  = tsfitting(data.site, data.lon, data.lat, data.decyr, data.E, data.SE, param_dict, 'E', timespan)
            eparam = erun.doFitting()

            # Vertical
            if dict_param['eqoffset_up'] == False: param_dict['eqlist']=[]
            if dict_param['eqpost_up'] == False: param_dict['eqpostlist']=[]
            urun  = tsfitting(data.site, data.lon, data.lat, data.decyr, data.U, data.SU, param_dict, 'U', timespan)
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
            if len(dict_output['eqpostdisp']) > 0:
                if len(dict_output['eqpost_tspan']) != 0:
                    output_postseismic_disp(nrun, erun, urun, dict_output['eqpost_tspan'], fid=fid_post_disp)
            if dict_output['eqpostts'] == True:
                if len(dict_output['eqpost_tspan']) != 0:
                    output_postseismic_ts(nrun, erun, urun, dict_output['eqpost_tspan'])
    if fid_velo      != None: fid_velo.close()
    if fid_eq        != None: fid_eq.close()
    if fid_brk       != None: fid_brk.close()
    if fid_post_disp != None: fid_post_disp.close()

if __name__ == '__main__':
    main()
