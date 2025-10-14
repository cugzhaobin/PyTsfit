import glob, sys
import numpy as np
from PyTsfit import eqcatalog, breakcatalog, eqPostList, correction, posData, tsfitting

eqfile     = './eq_rename.cors'
velfile    = './velomodel.vel.gmtvec'
velfile    = ''
offsetfile = ''
periodfile = ''
sitefile   = 'cmonoc.cmnc'
eq         = eqcatalog(eqfile)
bk         = breakcatalog(eqfile)
eqp        = eqPostList(eqfile, eq)
cor        = correction(velfile, offsetfile, periodfile)
plot_dict  = {'detrend':True, 'debreak': True, 'deeqoffset': True, 'depost': False, 'deseason':True}
mod_dict   = {'detrend':False, 'debreak': True, 'deeqoffset': True, 'depost': False, 'deseason':True}
sitelist   = np.genfromtxt(sitefile, dtype=str)

param_dict = {
              'constant'  : True,
              'linear'    : True,
              'eqlist'    : eq.eqlist,
              'eqpostlist': [],
              'brklist'   : bk.breaklist,
              'ANN'       : True,
              'SANN'      : True,
              'correct'   : cor}

tspan   = [2009, 2022]
fid     = open('cmnc_velo.gmtvec', 'a')
for i in range(len(sitelist)):
    posfiles = glob.glob('../pos/'+sitelist[i]+'*.pos')
    for j in range(len(posfiles)):
        posfile = posfiles[j]
        data = posData(posfile)

    # North
        param_dict['eqlist']     = eq.eqlist
        param_dict['eqpostlist'] = eqp.eqpostlist
        nrun  = tsfitting(data.site, data.lon, data.lat, data.decyr, data.N, data.SN, param_dict, 'N', 
                          tspan)
        nparam = nrun.doFitting()
        print(nparam)
    	
    	
        # East
        param_dict['eqlist']     = eq.eqlist
        param_dict['eqpostlist'] = eqp.eqpostlist
        erun  = tsfitting(data.site, data.lon, data.lat, data.decyr, data.E, data.SE, param_dict, 'E',
                          tspan)
        eparam = erun.doFitting()

        # Vertical
        param_dict['eqlist']     = []
        param_dict['eqpostlist'] = []
        urun  = tsfitting(data.site, data.lon, data.lat, data.decyr, data.U, data.SU, param_dict, 'U',
                          tspan)
        uparam = urun.doFitting()

#       plot_obs_mod(nrun, erun, urun, nparam, eparam, uparam, plot_dict, showFig=False)
        output_obs_mod(nrun, erun, urun, nparam, eparam, uparam, mod_dict)
#       output_velo(nrun, erun, urun, fid=fid)
#       output_eqoffset(nrun, erun, urun, fid=fid)

#       output_postseismic_disp(nrun, erun, urun, [2015.932, 2016.932], fid=fid)
