#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 19:48:51 2020

@author: zhao
"""

from PyTsfit import *
import glob, sys, time, os
import emcee, corner
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool

#os.environ["OMP_NUM_THREADS"] = "1"

def log_prior(theta, args):
    '''
    Calculate prior probability.
    
    Input:
        theta = variable
        args  = a list of input parameters
    Output:
        ln of prior probability.
    '''

    logic = [min(args[i])<theta[i]<max(args[i]) for i in range(len(args))]
    
    if sum(logic) == len(args):
        return 0.0
    else:
        return -np.inf


def log_likelihood(theta, run):
    '''
    Calculate likelihood.
    Mod by Zhao Bin, Jul. 30, 2019. Consider the co-variance matrix of InSAR data.
    
    Input:
        theta = variable
    Output:
        ln of likelihood.
    '''
  
    ifun = run.full_filter(run.t)
    obs  = run.obs
    mod  = ifun(run.t, *theta)
    res  = (obs-mod).reshape(len(obs),1)    
    cov  = run.sigma**2*np.eye(len(obs))
    icov = np.linalg.inv(cov)
#   print(res.T.dot(icov).dot(res))
    
    return -0.5*res.T.dot(icov).dot(res)

def log_posterior(theta, args, run):
    '''
    Calculate posterior probability based on prior distribution and likelihood distribution.
    
    Input:
        theta = variable
        leng  = min and max values of fault length
        wid   = min and max values of fault width
        dep   = min and max values of top depth of fault
        dip   = min and max values of fault dip
        strk  = min and max values of fault strike
        lat   = min and max values of fault location lat
        lon   = min and max values of fault location lon
    Output:
        ln of posterior probability.
    '''
    
    logic = [min(args[i])<theta[i]<max(args[i]) for i in range(len(args))]
    
    if sum(logic) == len(args):
        return log_prior(theta, args)+log_likelihood(theta, run)
    else:
        return -np.inf

eqfile     = './eq_rename.eq'
velfile    = './velomodel.vel.gmtvec'
velfile    = ''
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
              'ANN'       : False,
              'SANN'      : False,
              'correct'   : cor}


poslist = glob.glob('../pos/J041*.pos')
for posfile in poslist:
    data = posData(posfile)

    # North
    param_dict['eqlist']     = eq.eqlist
    param_dict['eqpostlist'] = eqp.eqpostlist
    nrun  = tsfitting(data.site, data.lon, data.lat, data.decyr, data.E, data.SE, param_dict, 'N', 
                 [2000, 2020])

    flag     = nrun.flag2
    ndim     = len(flag)
    popt     = []
    nwalkers = 20
    nburn    = 3000
    nsteps   = 5000

    for i in range(ndim):
        if flag[i] == 'CONSTANT':
            popt.append([-500, 500])
        if flag[i] == 'VELOCITY':
            popt.append([-100, 100])
        if flag[i] == 'EQOFFSET':
            popt.append([-5000, 5000])
        if flag[i] == 'EQDECAY':
            popt.append([-5000, 5000])
        if flag[i] == 'TAU':
            popt.append([5, 35])
        if flag[i] == 'BREAK':
            popt.append([-5000, 5000])
        if flag[i] == 'ANNUAL_SIN':
            popt.append([-20, 20])
        if flag[i] == 'ANNUAL_COS':
            popt.append([-20, 20])
        if flag[i] == 'SANNUAL_SIN':
            popt.append([-20, 20])
        if flag[i] == 'SANNUAL_COS':
            popt.append([-20, 20])
    popt = np.array(popt)
    print(ndim, len(popt))
    
    starting_guess = np.random.random((nwalkers, ndim))
    for i in range(ndim):
        starting_guess[:,i] = np.random.uniform(min(popt[i]), max(popt[i]), nwalkers)
#   with Pool() as pool:
#       sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, pool=pool, args=[popt, nrun])
#       start = time.time()
#       sampler.run_mcmc(starting_guess, nsteps)
#       end   = time.time()
#       emcee_trace = sampler.chain[:,nburn:,:].reshape(-1, ndim)
#       print(end-start)
#       fig = corner.corner(emcee_trace, show_titles=True)
#       fig.savefig("posterior_distribution.png")
#       np.save('mcmc_trace', emcee_trace)
#   
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, threads=1, args=[popt, nrun])
    start = time.time()
    sampler.run_mcmc(starting_guess, nsteps)
    end   = time.time()
    emcee_trace = sampler.chain[:,nburn:,:].reshape(-1, ndim)
    print(end-start)
    fig = corner.corner(emcee_trace, show_titles=True)
    fig.savefig("posterior_distribution.png")
    np.save('mcmc_trace', emcee_trace)
    
