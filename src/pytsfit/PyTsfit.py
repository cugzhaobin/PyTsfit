#!/usr/bin/env python
# ----------------------------------------------------------
# PyTsfit: fitting GPS coordinate time series or baseline time series
#
# Zhao Bin, Institute of Seismologym CEA.
# Nov 12, 2018
#
# This file is now a thin compatibility shim. The implementation lives in
# data.py / models.py / tsfitting.py / output.py; this module re-exports the
# same names so ``from PyTsfit import *`` and ``from pytsfit.PyTsfit import *``
# keep working after the split.
# ----------------------------------------------------------
import os, sys, logging
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

logging.basicConfig(level=logging.INFO,
    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
    datefmt="%d-%m-%Y %H:%M:%S")

from pytsfit import *  # noqa: F401,F403  (re-export the public API)

# Re-export the names the old module exposed at top level (beyond its own
# classes/functions), in case downstream ``import *`` relies on them.
__all__ = (['gt', 'gpstime', 'neuData', 'posData',
            'earthquake', 'eqcatalog', 'eqPost', 'eqPostList',
            'offset', 'breakcatalog', 'correction', 'tsfitting',
            'output_velo', 'output_postseismic_velo', 'output_postseismic_disp',
            'output_eqoffset', 'output_break', 'output_postseismic_ts',
            'output_period', 'output_summary', 'plot_obs_mod',
            'output_obs_mod', 'output_param']
           + ['np', 'plt', 'curve_fit', 'os', 'sys', 'logging'])


if __name__ == '__main__':
    import time, getpass
    print("# Created by %s on %s" %(getpass.getuser(), time.asctime()))
    text = "from PyTsfit import *\nimport glob, os\n\neqfile     = ''\nvelfile    = ''\noffsetfile = ''\nperiodfile = ''\n"
    print("%s" %(text))
    text = "eq         = eqcatalg(eqfile)\nbk         = breakcatalog(eqfile)\neqp        = eqPostList(eqfile, eq)"
    print("%s" %(text))
    text = "cor        = correction(velfile, offsetfile, periodfile\n"
    print("%s" %(text))
    text = "plot_dict  = {'detrend':True, 'debreak':True, 'deeqoffset':False, 'depost':False, 'deseason':False}"
    print("%s" %(text))
    text = "mod_dict   = {'detrend':True, 'debreak':True, 'deeqoffset':False, 'depost':False, 'deseason':False}"
    print("%s" %(text))
    text = "param_dict = {\n             'constant'   : True,\n             'linear'     : True,"
    print("%s" %(text))
    text = "             'eqlist'     : eq.eqlist,\n             'eqpostlist' : [],\n             'brklist'    : bk.breaklist"
    print("%s" %(text))
    text = "             'ANN'        : False,\n             'SANN'       : False\n             'correct'    : cor}\n"
    print("%s" %(text))
    print("poslist  = glob.glob('*.pos')")
    print("fid      = open('velo.gmtvec', 'a')")
    print("for posfile in poslist:")
    print("    data = posData(posfile)")
    print("    param_dict['eqpostlist'] = eqp.eqpostlist")
    print("    nrun = tsfitting(data.site, data.lon, data.lat, data.decyr, data.N, data.SN, param_dict, 'N', [2000,2020])")
    print("    nparam = nrun.doFitting()\n")
    print("    param_dict['eqpostlist'] = eqp.eqpostlist")
    print("    erun = tsfitting(data.site, data.lon, data.lat, data.decyr, data.E, data.SE, param_dict, 'E', [2000,2020])")
    print("    eparam = erun.doFitting()\n")
    print("    param_dict['eqpostlist'] = []")
    print("    urun = tsfitting(data.site, data.lon, data.lat, data.decyr, data.U, data.SU, param_dict, 'U', [2000,2020])")
    print("    uparam = erun.doFitting()\n")
    print("    plot_obs_mod(nrun, erun, urun, nparam, eparam, uparam, plot_dict)")
    print("    output_postseismic_disp(nrun, erun, urun, [2015.315, 2017.315], fid=fid)")
    print("    output_postseismic_velo(nrun, erun, urun, [2015.315, 2017.315], fid=fid)")
    print("    output_velo(nrun, erun, urun, fid=fid)")
    print("    output_eqoffset(nrun, erun, urun)")
    print("    output_break(nrun, erun, urun)")
