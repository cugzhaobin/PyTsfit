#from PyTsfit.tsfitting import tsfitting
from PyTsfit import *
import glob
import sys

eqfile     = './eq_rename.eq'
velfile    = ''
velfile    = './velomodel.vel.gmtvec'
offsetfile = ''
periodfile = ''
eq         = eqcatalog(eqfile)
bk         = breakcatalog(eqfile)
eqp        = eqPostList(eqfile, eq)
cor        = correction(velfile, offsetfile, periodfile)
plot_dict  = {'detrend':True, 'debreak': True, 'deeqoffset': False, 'depost': False, 'deseason':False}
mod_dict   = {'detrend':False, 'debreak': False, 'deeqoffset': False, 'depost': True, 'deseason':False}


param_dict = {
              'constant'  : True,
              'linear'    : True,
              'eqlist'    : eq.eqlist,
              'eqpostlist': [],
              'brklist'   : bk.breaklist,
              'ANN'       : True,
              'SANN'      : True,
              'correct'   : cor}


poslist = glob.glob('*.pos')
fid     = open('gorkha_2yr.gmtvec', 'a')
for posfile in poslist:
    print(posfile)
    data = posData(posfile)

    # North
    param_dict['eqlist']     = eq.eqlist
    param_dict['eqpostlist'] = eqp.eqpostlist
    nrun  = tsfitting(data.site, data.lon, data.lat, data.decyr, data.N, data.SN, param_dict, 'N', 
                 [2000, 2020])
    nparam = nrun.doFitting()


    # East
    param_dict['eqlist']     = eq.eqlist
    param_dict['eqpostlist'] = eqp.eqpostlist
    erun  = tsfitting(data.site, data.lon, data.lat, data.decyr, data.E, data.SE, param_dict, 'E',
                 [2000, 2020])
    eparam = erun.doFitting()
    print(eparam)

    # Vertical
    param_dict['eqlist']     = []
    param_dict['eqpostlist'] = []
    urun  = tsfitting(data.site, data.lon, data.lat, data.decyr, data.U, data.SU, param_dict, 'U',
                 [2000, 2020])
    uparam = urun.doFitting()

    plot_obs_mod(nrun, erun, urun, nparam, eparam, uparam, plot_dict)
#   output_postseismic_disp(nrun, erun, urun, mod_dict, [2015.315, 2017.315], fid=fid)
