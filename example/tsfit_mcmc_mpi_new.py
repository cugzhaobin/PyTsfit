#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 19:48:51 2020

@author: zhao
"""

from PyTsfit import *
import glob, sys, time
import emcee, corner, argparse
import numpy as np
from schwimmbad import MPIPool


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
    Mod by Zhao Bin, Jul. 30, 2019. Compute likelihood
    
    Input:
        theta = variable
    Output:
        ln of likelihood.
    '''
  
    ifun = run.full_filter(run.t)
    obs  = run.obs
    mod  = ifun(run.t, *theta)
#   res  = (obs-mod).reshape(len(obs),1)    
#   cov  = run.sigma**2*np.eye(len(obs))
#   icov = np.linalg.inv(cov)
#   return -0.5*res.T.dot(icov).dot(res)
    res  = (obs-mod)/run.sigma
    return -0.5*res.T.dot(res)

def log_posterior(theta, args, run):
    '''
    Calculate posterior probability based on prior distribution and likelihood distribution.
    
    Input:
        theta = variable
        args  = a list of input parameters
    Output:
        ln of posterior probability.
    '''
    
    logic = [min(args[i])<theta[i]<max(args[i]) for i in range(len(args))]
    
    if sum(logic) == len(args):
        return log_prior(theta, args)+log_likelihood(theta, run)
    else:
        return -np.inf

def set_bound(flag, cor, site, component='E'):
    '''
    '''
    popt = []
    ndim = len(flag)
    for i in range(ndim):
        if flag[i] == 'CONSTANT':
            popt.append([-500, 500])
        if flag[i] == 'VELOCITY':
            if cor.correct_velo == True:
                if len(np.where(cor.velsite == site)[0]) > 0:
                    if component == 'E':
                        ve = cor.velsite[cor.velsite==site][0,2]
                        popt.append([ve-1.0, ve+1.0])
                    elif component == 'N':
                        vn = cor.velsite[cor.velsite==site][0,3]
                        popt.append([ve-1.0, ve+1.0])
                    else:
                        popt.append([-100, 100])
                else:
                    popt.append([-100, 100])
            else:
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
    return popt

def main(args):
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
    timespan   = [2011, 2025]
    nburns     = args.nburns
    nsteps     = args.nsteps
    
    param_dict = {
                  'constant'  : True,
                  'linear'    : True,
                  'eqlist'    : eq.eqlist,
                  'eqpostlist': [],
                  'brklist'   : bk.breaklist,
                  'ANN'       : False,
                  'SANN'      : False,
                  'correct'   : cor}


    poslist = glob.glob('./2837*.pos')
    for posfile in poslist:
        data = posData(posfile)

        # East component
        param_dict['eqlist'] = eq.eqlist
        erun  = tsfitting(data.site, data.lon, data.lat, data.decyr, data.E, data.SE, param_dict, 'E', timespan)
        flag  = erun.flag2
        ndim  = len(flag)
        nwalkers = 2*ndim
        popt  = set_bound(flag, cor, data.site, component = 'E')
        print(flag)
    
        starting_guess = np.random.random((nwalkers, ndim))
        for i in range(ndim):
            starting_guess[:,i] = np.random.uniform(min(popt[i]), max(popt[i]), nwalkers)
    
        with MPIPool() as pool:
            if not pool.is_master():
                pool.wait()
                sys.exit(0)
            sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, pool=pool, args=[popt, erun])
            sampler.run_mcmc(starting_guess, nsteps, progress=True)
            chain   = sampler.get_chain()
            np.savez('chain', chain)
            emcee_trace = chain[nburns:,:,:].reshape(-1, ndim)
            fig         = corner.corner(emcee_trace, show_titles=True)
            fig.savefig("{}_{}_posterior.png".format(data.site, "E"))


        # North component
        param_dict['eqlist'] = eq.eqlist
        nrun  = tsfitting(data.site, data.lon, data.lat, data.decyr, data.N, data.SN, param_dict, 'N', timespan)
        flag  = nrun.flag2
        ndim  = len(flag)
        nwalkers = 2*ndim
        popt  = set_bound(flag, cor, data.site, component = 'N')
        print(flag)
    
        starting_guess = np.random.random((nwalkers, ndim))
        for i in range(ndim):
            starting_guess[:,i] = np.random.uniform(min(popt[i]), max(popt[i]), nwalkers)
    
        with MPIPool() as pool2:
            if not pool2.is_master():
                pool2.wait()
                sys.exit(0)
            sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, pool=pool2, args=[popt, nrun])
            sampler.run_mcmc(starting_guess, nsteps, progress=True)
            chain   = sampler.get_chain()
            np.savez('{}_{}_chain'.format(data.site, "N"), chain)
            emcee_trace = chain[nburns:,:,:].reshape(-1, ndim)
            fig         = corner.corner(emcee_trace, show_titles=True)
            fig.savefig("{}_{}_posterior.png".format(data.site, "N"))
        # Up component
#       urun  = tsfitting(data.site, data.lon, data.lat, data.decyr, data.U, data.SU, param_dict, 'U', timespan)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Fit time series using MCMC method")
    parser.add_argument('--nsteps', type=str, required=False, default=2000)
    parser.add_argument('--nburns',  type=str, required=False, default=1000)
    parser.add_argument('--nwalkers', type=str, required=False, default=10)
    args = parser.parse_args()
    main(args)
