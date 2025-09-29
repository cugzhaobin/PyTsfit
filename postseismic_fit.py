#!/usr/bin/env python
#from PyTsfit.tsfitting import tsfitting
from PyTsfit import *
import glob
import sys
import numpy as np

eqfile     = './eq_rename.eq'
velfile    = ''
velfile    = './velomodel.vel.gmtvec'
offsetfile = ''
periodfile = './season_new.dat'
sitefile   = 'sitelist'
eq         = eqcatalog(eqfile)
bk         = breakcatalog(eqfile)
eqp        = eqPostList(eqfile, eq)
cor        = correction(velfile, offsetfile, periodfile)
plot_dict  = {'detrend':True, 'debreak': True, 'deeqoffset': True, 'depost': False, 'deseason':False}
mod_dict   = {'detrend':False, 'debreak': False, 'deeqoffset': False, 'depost': True, 'deseason':False}
sitelist   = np.genfromtxt(sitefile, dtype='S')

param_dict = {
              'constant'  : True,
              'linear'    : True,
              'eqlist'    : eq.eqlist,
              'eqpostlist': [],
              'brklist'   : bk.breaklist,
              'ANN'       : False,
              'SANN'      : False,
              'correct'   : cor}


fid     = open('velocity.gmtvec', 'a')
for i in range(len(sitelist)):
    posfiles = glob.glob('./'+sitelist[i].decode()+'*.pos')
    for j in range(len(posfiles)):
	    posfile = posfiles[j]
	    print(posfile)
	    data = posData(posfile)
	
	    # North
	    param_dict['eqlist']     = eq.eqlist
	    param_dict['eqpostlist'] = eqp.eqpostlist
	    nrun  = tsfitting(data.site, data.lon, data.lat, data.decyr, data.N, data.SN, param_dict, 'N', 
	                 [1998, 2021])
	    nparam = nrun.doFitting()
	    print(nparam)
	
	
	    # East
	    param_dict['eqlist']     = eq.eqlist
	    param_dict['eqpostlist'] = eqp.eqpostlist
	    erun  = tsfitting(data.site, data.lon, data.lat, data.decyr, data.E, data.SE, param_dict, 'E',
	                 [1998, 2021])
	    eparam = erun.doFitting()
	
	    # Vertical
	    param_dict['eqlist']     = []
	    param_dict['eqpostlist'] = []
	    urun  = tsfitting(data.site, data.lon, data.lat, data.decyr, data.U, data.SU, param_dict, 'U',
	                 [1998, 2021])
	    uparam = urun.doFitting()
	    output_postseismic_disp(nrun, erun, urun, [2015.315, 2019.315], fid=fid)
	    plot_obs_mod(nrun, erun, urun, nparam, eparam, uparam, plot_dict)
