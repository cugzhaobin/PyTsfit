#!/usr/bin/env python
# ----------------------------------------------------------
# The time series fitting engine.
#
# Extracted from PyTsfit.py during the module split. No behaviour change.
# ----------------------------------------------------------
import os, sys, logging
import numpy as np
import matplotlib.pyplot as plt
from . import geotools as gt
from scipy.optimize import curve_fit

class tsfitting:
    '''
    tsfitting is a class representing fitting a signle time series

    Mod by Zhao Bin, Jan. 9, 2019. Fix bug sys.exit() --> return
    '''

    # number of parameter
    nparam = 0

    # decimal year
    t      = np.array([])

    # observation
    obs    = np.array([])

    # sigma
    sigma  = np.array([])

    # flag
    flag   = []

    # earthquake list
    ieqlist  = []

    # post earthquake list
    ieqpostlist = []

    # break list
    ibrklist = []

    # parmeters
    param = []

    def __init__(self, site, lon, lat, t, obs, sigma, param_dict, component='N', time_range=None):
        '''
        Constructor.

        Input:
            site       = 4 char site ID
            lon        = longitude
            lat        = latitude
            t          = time list
            obs        = observation list
            sigma      = uncertainty of observation
            param_dict = {}
            component  = "N/E/U"
            time_range = time range for fitting
        '''
        self.nparam   = 0
        self.t        = []
        self.obs      = []
        self.sigma    = []
        self.site     = ""
        self.midt     = None
        self.flag     = []
        self.flag2    = []
        self.param    = []
        self.ieqlist  = []
        self.ieqpostlist = []
        self.ibrklist = []
        self.component= component

        if time_range is None:
            time_range = [-np.inf, np.inf]
        idx  = np.where(np.logical_and(t>time_range[0], t<time_range[1]))[0]
        if len(idx) < 2: return
        self.lon      = lon
        self.lat      = lat
        self.t        = t[idx]
        self.obs      = obs[idx]
        self.sigma    = sigma[idx]
        self.site     = site
        self.midt     = (self.t[0]+self.t[-1])/2
        mint          = min(self.t)
        maxt          = max(self.t)


        # Constant term
        if param_dict['constant'] == True:
            self.nparam = self.nparam + 1
            self.flag.append('CONSTANT')
            self.flag2.append('CONSTANT')
#           print('constant: %d' %(self.nparam))

        # Linear term
        if param_dict['linear'] == True:
            self.nparam = self.nparam + 1
            self.flag.append('VELOCITY')
            self.flag2.append('VELOCITY')
#           print('velocity: %d' %(self.nparam))

        # Coseismic term
        for eq in param_dict['eqlist']:
            origin = eq.location
            xy   = gt.llh2localxy([lat, lon], origin)
            dist = np.sqrt(xy[0]**2 + xy[1]**2)

            if dist < eq.distance and eq.decyr > mint and eq.decyr < maxt:
                self.ieqlist.append(eq)
                self.nparam = self.nparam + 1
                self.flag.append('EQOFFSET')
                self.flag2.append('EQOFFSET')
#               print('eqoffset: %d, %s' %(self.nparam, eq.code))

        # Postseismic term
        for eqpost in param_dict['eqpostlist']:
            origin = eqpost.eq.location
            xy     = gt.llh2localxy([lat, lon], origin)
            dist   = np.sqrt(xy[0]**2 + xy[1]**2)
            if dist < eqpost.eq.distance and eqpost.eq.decyr < maxt:
                self.ieqpostlist.append(eqpost)
                self.nparam = self.nparam + 2
                self.flag.append("EQDECAY")
                self.flag2.append("EQDECAY")
                self.flag2.append("TAU")
#               print('eqpost: %d' %(self.nparam))

        # Break term
        for brk in param_dict['brklist']:
            if brk.site == site and brk.decyr > mint and brk.decyr < maxt:
                self.ibrklist.append(brk)
                self.nparam = self.nparam + 1
                self.flag.append('BREAK')
                self.flag2.append('BREAK')
#               print('break: %d' %(self.nparam))

        # Period term - annual
        if param_dict['ANN'] == True:
            self.nparam = self.nparam + 2
            self.flag.append('ANNUAL')
            self.flag2.append('ANNUAL_SIN')
            self.flag2.append('ANNUAL_COS')
#           print('A period: %d' %(self.nparam))

        # Period term - semi-annual
        if param_dict['SANN'] == True:
            self.nparam = self.nparam + 2
            self.flag.append('SANNUAL')
            self.flag2.append('SANNUAL_SIN')
            self.flag2.append('SANNUAL_COS')
#           print('S period: %d' %(self.nparam))

        # correction
        if 'correct' in param_dict.keys():
            self.correct = param_dict['correct']

        self.flag  = np.array(self.flag)
        self.flag2 = np.array(self.flag2)
        param = ' '.join(self.flag2)
        logging.info('Parameters {}'.format(param))

    def setBoundAndInit(self):
        '''
        Set the bounds and init values for curve_fit function

        Mod by Zhao Bin, Jan 21, 2019. Adding constraints on EQOFFSET and BREAK
        Mod by Zhao Bin, Feb 17, 2019. Fix bug for period constraints
        Mod by Zhao Bin, Jun 10, 2020. Fix bug for EQOFFSET
        '''
        lb    = [-np.inf]*self.nparam
        ub    =  [np.inf]*self.nparam
        pinit = [10]*self.nparam
        eqid  = 0
        brkid = 0

        i, j, k = 0, 0, 0
        for item in self.flag:
            if item == 'CONSTANT':
                i = i+1
            if item == 'VELOCITY':
                if self.correct.correct_velo:
                    if self.component == 'N':
                        if len(np.where(self.correct.velsite==self.site)[0])>0:
                            pinit[i]  = self.correct.veldata[self.correct.velsite==self.site][0,3]
                            lb[i]     = pinit[i]-1e-1
                            ub[i]     = pinit[i]+1e-1
                        else:
                            logging.warning('No prior velocity for site {}'.format(self.site))
                    elif self.component == 'E':
                        if len(np.where(self.correct.velsite==self.site)[0])>0:
                            pinit[i]  = self.correct.veldata[self.correct.velsite==self.site][0,2]
                            lb[i]     = pinit[i]-1e-1
                            ub[i]     = pinit[i]+1e-1
                        else:
                            logging.warning('No prior velocity for site {}'.format(self.site))
                    elif self.component == 'U':
                        if self.correct.veldata.shape[1] <= 8:
                            lb[i]     = -np.inf
                            ub[i]     = np.inf
                        else:
                            if len(np.where(self.correct.velsite==self.site)[0])>0:
                                pinit[i]  = self.correct.veldata[self.correct.velsite==self.site][0,8]
                    else:
                        pinit[i]  = 0.0
                i = i+1
            if item == 'EQOFFSET':
                if self.correct.correct_offset:
                    for m in range(eqid, len(self.ieqlist)):
                        eq   = self.ieqlist[m]
                        eqid = eqid+1
                        if len(np.where(self.correct.offsetsite==self.site)[0])>0 and\
                                len(np.where(abs(self.correct.offsetyear-eq.decyr)<0.01)[0]):
                            indx1 = set(np.where(self.correct.offsetsite==self.site)[0])
                            indx2 = set(np.where(abs(self.correct.offsetyear-eq.decyr)<0.01)[0])
                            indx  = list(indx1 & indx2)
                            if self.component == 'N' and len(indx)>0:
                                pinit[i] = self.correct.offsetdata[indx][0,1]
                                lb[i]    = pinit[i]-1e-1
                                ub[i]    = pinit[i]+1e-1
                                i        = i+1
                                j        = j+1
                                continue
                            elif self.component == 'E' and len(indx)>0:
                                pinit[i] = self.correct.offsetdata[indx][0,0]
                                lb[i]    = pinit[i]-1e-1
                                ub[i]    = pinit[i]+1e-1
                                i        = i+1
                                j        = j+1
                                continue
                            elif self.component == 'U' and len(indx)>0:
                                pinit[i] = self.correct.offsetdata[indx][0,2]
                                lb[i]    = pinit[i]-1e-1
                                ub[i]    = pinit[i]+1e-1
                                i        = i+1
                                j        = j+1
                                continue
                        else:
                            i = i+1
                            j = j+1
                else:
                    i = i+1
                    j = j+1
            if item == 'EQDECAY':
                lb[i+1]    = self.ieqpostlist[k].mintau
                ub[i+1]    = self.ieqpostlist[k].maxtau
                pinit[i+1] = lb[i+1]
                i = i+2
                j = k+1
            if item == 'BREAK':
                if self.correct.correct_offset:
                    for m in range(brkid, len(self.ibrklist)):
                        brk   = self.ibrklist[m]
                        brkid = brkid+1
                        if len(np.where(self.correct.offsetsite==self.site)[0])>0 and\
                                len(np.where(abs(self.correct.offsetyear-brk.decyr)<0.01)[0]):
                            indx1 = set(np.where(self.correct.offsetsite==self.site)[0])
                            indx2 = set(np.where(abs(self.correct.offsetyear-brk.decyr)<0.01)[0])
                            indx  = list(indx1 & indx2)
                            if self.component == 'N' and len(indx)>0:
                                pinit[i] = self.correct.offsetdata[indx][0,1]
                                lb[i]    = pinit[i]-1e-1
                                ub[i]    = pinit[i]+1e-1
                                break
                            elif self.component == 'E' and len(indx)>0:
                                pinit[i] = self.correct.offsetdata[indx][0,0]
                                lb[i]    = pinit[i]-1e-1
                                ub[i]    = pinit[i]+1e-1
                                break
                            elif self.component == 'U' and len(indx)>0:
                                pinit[i] = self.correct.offsetdata[indx][0,2]
                                lb[i]    = pinit[i]-1e-1
                                ub[i]    = pinit[i]+1e-1
                                break
                i = i+1

            if item == 'ANNUAL':
                if self.correct.correct_period:
                    idx = np.where(self.correct.periodsite==self.site)[0]
                    if self.component == 'N':
                        if len(idx) == 0:
                            pinit[i], pinit[i+1]  = 0.0, 0.0
                        else:
                            pinit[i]   = self.correct.perioddata[idx][0,4]
                            pinit[i+1] = self.correct.perioddata[idx][0,5]
                    elif self.component == 'E':
                        if len(idx)==0:
                            pinit[i], pinit[i+1]  = 0.0, 0.0
                        else:
                            pinit[i]   = self.correct.perioddata[idx][0,0]
                            pinit[i+1] = self.correct.perioddata[idx][0,1]
                    elif self.component == 'U':
                        if len(idx)==0:
                            pinit[i], pinit[i+1]  = 0.0, 0.0
                            logging.warning('No proior annual term for site {}'.format(self.site))
                        else:
                            pinit[i]   = self.correct.perioddata[idx][0,8]
                            pinit[i+1] = self.correct.perioddata[idx][0,9]
                    lb[i]      = pinit[i]-1e-1
                    ub[i]      = pinit[i]+1e-1
                    lb[i+1]    = pinit[i+1]-1e-1
                    ub[i+1]    = pinit[i+1]+1e-1
                i = i+2
            if item == 'SANNUAL':
                if self.correct.correct_period:
                    idx = np.where(self.correct.periodsite==self.site)[0]
                    if self.component == 'N':
                        if len(idx)==0:
                            pinit[i], pinit[i+1]  = 0.0, 0.0
                        else:
                            pinit[i]   = self.correct.perioddata[idx][0,6]
                            pinit[i+1] = self.correct.perioddata[idx][0,7]
                    elif self.component == 'E':
                        if len(idx)==0:
                            pinit[i], pinit[i+1]  = 0.0, 0.0
                        else:
                            pinit[i]   = self.correct.perioddata[idx][0,2]
                            pinit[i+1] = self.correct.perioddata[idx][0,3]
                    elif self.component == 'U':
                        if len(idx)==0:
                            pinit[i], pinit[i+1]  = 0.0, 0.0
                            logging.info('No proior semi-annual term for site {}'.format(self.site))
                        else:
                            pinit[i]   = self.correct.perioddata[idx][0,10]
                            pinit[i+1] = self.correct.perioddata[idx][0,11]
                    lb[i]      = pinit[i]-1e-1
                    ub[i]      = pinit[i]+1e-1
                    lb[i+1]    = pinit[i+1]-1e-1
                    ub[i+1]    = pinit[i+1]+1e-1
                i = i+2

        return lb, ub, pinit

    def doFitting(self):
        '''
        Do fitting using the ifun
        Mod by Zhao Bin, Jan. 10, 2019. Fix bug when no data to be fitted and compute WRMS

        Output:
            popt  = estiated parameters
        '''

        if len(self.flag) == 0: return np.empty(0)

        # function
        ifun          = self.full_filter(self.t)

        # set lower and upper bounds and initial parameters
        lb, ub, pinit = self.setBoundAndInit()

        # fit curve
        try:
            popt, pcov = curve_fit(ifun,
                               self.t,
                               self.obs,
                               sigma=self.sigma,
                               p0=pinit,
                               bounds=[lb, ub])
            self.param = popt
            self.cov   = pcov
            self.ifun  = ifun
            self.res   = self.obs-ifun(self.t, *popt)
            self.wrms  = np.sqrt(sum((self.res/self.sigma)**2)/sum(1.0/self.sigma**2))
            return popt
        except:
            return np.array([])


    def full_filter(self, t):
        '''
        Construct a function depending on flag

        Input:
            t          = a list of time

        Output:
            pos_filter = function
        '''

        def pos_filter(t, *p):
            '''
            Inner function

            Input:
                t  = time list
                p  = arguments
            '''
            y = 0
            i = 0
            j = 0
            k = 0
            l = 0
            for item in self.flag:
                if item == 'CONSTANT':
                    y  = y + p[i]
                    i  = i + 1
                if item == 'VELOCITY':
                    y  = y + p[i] * (t-self.midt)
                    i  = i + 1
                if item == 'EQOFFSET':
                    t0 = self.ieqlist[j].decyr
                    y  = y + p[i]*np.heaviside(t - t0, 0)
                    i  = i + 1
                    j  = j + 1
                if item == 'EQDECAY':
                    t0 = self.ieqpostlist[k].eq.decyr
                    if self.ieqpostlist[k].method == 'LOG':
                        temp = np.zeros(len(t))
                        idx  = np.where(t>t0)
                        temp[idx] = p[i]*np.log(1+(t[idx]-t0)*365.25/p[i+1])
                        y = y + temp
                        i = i + 2
                    if self.ieqpostlist[k].method == 'EXP':
                        y = y + np.heaviside(t - t0, 0)*\
                            p[i]*(1-np.exp(-(t-t0)*365.25/p[i+1]))
                        i = i + 2
                    k  = k + 1
                if item == 'BREAK':
                    t0 = self.ibrklist[l].decyr
                    y  = y + p[i]*np.heaviside(t - t0, 0)
                    i  = i + 1
                    l  = l + 1
                if item == 'ANNUAL':
                    y  = y + p[i]  * np.sin(2*np.pi*t) +\
                             p[i+1]* np.cos(2*np.pi*t)
                    i  = i + 2
                if item == 'SANNUAL':
                    y  = y + p[i]  * np.sin(4*np.pi*t) +\
                             p[i+1]* np.cos(4*np.pi*t)
                    i  = i + 2
            return y
        return pos_filter

    def plot_obs_mod(self):
        '''
        Plot time series of observation and models
        '''
        mt  = np.arange(min(self.t), max(self.t), 1/365.25)
        model = self.ifun(mt, *self.param)
        plt.plot(self.t, self.obs, 'ro', ms=2)
        plt.plot(self.t, model, color='b')

        plt.title(self.site)
        plt.xlabel('Time (year)')
        plt.ylabel('Displacement (mm)')
        plt.show()

    def get_mod(self, time_span=None):
        '''
        Output modeled time series

        Mod by Zhao Bin, Jan 31, 2019. Adding parameter time_span
        '''
        if time_span is not None and len(time_span) == 2:
            mt = np.arange(min(time_span), max(time_span), 1/365.25)
        else:
            mt = np.arange(min(self.t), max(self.t), 1/365.25)
        mfun = self.full_filter(mt)
        m    = mfun(mt, *self.param)
        return mt, m

    def get_correct(self, mod_dict, time_span=None, eqcode=None):
        '''
        Retrieve modeled time series according to mod_dict at observed and modeled time list

        Mod by Zhao Bin, Jan 31, 2019. Adding parameter time_span and fix a bug

        Input:
            mod_dict  = {}
            time_span = []
            eqcode    = earthquake code with 2char
        Output:
            obs_correct = modeled time series at observed epochs
            mod_correct = modeled time series at modeled epochs
        '''

        param = self.param.copy()
        if mod_dict['detrend'] == False:
            idx        = np.where(self.flag2 == 'VELOCITY')[0]
            param[idx] = 0.0
        if mod_dict['deeqoffset'] == False:
            idx        = np.where(self.flag2 == 'EQOFFSET')[0]
            param[idx] = 0.0
        if mod_dict['depost'] == False:
            idx        = np.where(self.flag2 == 'EQDECAY')[0]
            param[idx] = 0.0
        if mod_dict['debreak'] == False:
            idx        = np.where(self.flag2 == 'BREAK')[0]
            param[idx] = 0.0
        if mod_dict['deseason'] == False:
            idx        = np.where(np.logical_or(self.flag2=='ANNUAL_SIN', self.flag2=='ANNUAL_COS'))[0]
            param[idx] = 0.0
            idx        = np.where(np.logical_or(self.flag2=='SANNUAL_SIN', self.flag2=='SANNUAL_COS'))[0]
            param[idx] = 0.0

        # Added by Zhao Bin, Mar. 11, 2022
        for i, eqpost in enumerate(self.ieqpostlist):
            if eqpost.eq.code != eqcode and eqcode != None:
                idx    = np.where(self.flag2 == 'EQDECAY')[0]
                param[idx[i]] = 0.0
                logging.info('set postseismic amplitude of eq {} to zero'.format(eqpost.eq.code))

        obs_correct = np.array([])
        mod_correct = np.array([])
        if time_span is not None and len(time_span) == 2:
            idx    = np.where(np.logical_and(self.t>time_span[0], self.t<time_span[1]))[0]
            if len(idx) > 0:
                obs_correct = self.ifun(self.t[idx], *param)
        else:
            obs_correct = self.ifun(self.t, *param)
        if time_span is not None and len(time_span) == 2:
            mt = np.arange(min(time_span), max(time_span), 1/365.25)
        else:
            mt = np.arange(min(self.t), max(self.t), 1/365.25)
        mfun        = self.full_filter(mt)
        mod_correct = mfun(mt, *param)

        return obs_correct, mod_correct
