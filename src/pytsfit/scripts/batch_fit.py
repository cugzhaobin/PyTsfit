#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Batch fitting of PBO .pos time series.

Reads a site list and fits N/E/U components of every ``*.pos`` file found
under ``../pos/``, then writes the observed/modelled series and the secular
velocities. Update the file paths and parameter dictionary for your own data.

Written for the post-split package layout: names are imported directly from
the owning modules (data / models / tsfitting / output).
'''
import glob
import numpy as np
from pytsfit.data import posData
from pytsfit.models import eqcatalog, breakcatalog, eqPostList, correction
from pytsfit.tsfitting import tsfitting
from pytsfit.output import plot_obs_mod, output_obs_mod, output_velo

eqfile     = './eq_rename.cors'
velfile    = './velomodel.vel.gmtvec'
offsetfile = ''
periodfile = ''
sitefile   = 'cmonoc.cmnc'
eq         = eqcatalog(eqfile)
bk         = breakcatalog(eqfile)
eqp        = eqPostList(eqfile, eq)
cor        = correction(velfile, offsetfile, periodfile)

plot_dict  = {'detrend': True, 'debreak': True, 'deeqoffset': True,
              'depost': False, 'deseason': True,
              'figformat': 'jpg', 'showfig': False}
mod_dict   = {'detrend': False, 'debreak': True, 'deeqoffset': True,
              'depost': False, 'deseason': True}

# Optional fitting / quality-control options. Leave as None for the defaults,
# or pass a dict, e.g. {'sigma_scale': 'realistic', 'outlier': True}.
fit_opts   = None

sitelist   = np.genfromtxt(sitefile, dtype=str)
if sitelist.size == 1:
    sitelist = np.array([str(sitelist)])

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
    posfiles = glob.glob('../pos/' + sitelist[i] + '*.pos')
    for j in range(len(posfiles)):
        posfile = posfiles[j]
        data    = posData(posfile)

        # North
        param_dict['eqlist']     = eq.eqlist
        param_dict['eqpostlist'] = eqp.eqpostlist
        nrun   = tsfitting(data.site, data.lon, data.lat, data.decyr, data.N, data.SN,
                           param_dict, 'N', tspan, fit_opts=fit_opts)
        nparam = nrun.doFitting()
        print(nparam)

        # East
        param_dict['eqlist']     = eq.eqlist
        param_dict['eqpostlist'] = eqp.eqpostlist
        erun   = tsfitting(data.site, data.lon, data.lat, data.decyr, data.E, data.SE,
                           param_dict, 'E', tspan, fit_opts=fit_opts)
        eparam = erun.doFitting()

        # Vertical
        param_dict['eqlist']     = []
        param_dict['eqpostlist'] = []
        urun   = tsfitting(data.site, data.lon, data.lat, data.decyr, data.U, data.SU,
                           param_dict, 'U', tspan, fit_opts=fit_opts)
        uparam = urun.doFitting()

        plot_obs_mod(nrun, erun, urun, nparam, eparam, uparam, plot_dict)
        output_obs_mod(nrun, erun, urun, nparam, eparam, uparam, mod_dict)
        output_velo(nrun, erun, urun, fid=fid)

fid.close()
