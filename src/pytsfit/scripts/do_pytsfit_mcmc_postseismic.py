#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 19:48:51 2020

@author: zhao
"""

from PyTsfit import posData, tsfitting
import glob, sys, time
import emcee, corner, argparse
import numpy as np
from scipy.linalg import norm
import matplotlib.pyplot as plt


def postsesimic_trace(flag, trace, eqtime, timespan, method='LOG'):
    '''
    '''
    flag = list(flag)
    if method == 'LOG':
        tau = trace[:,flag.index('TAU')]
        amp = trace[:,flag.index('EQDECAY')]
        post = amp*np.log(1+((timespan[1]-eqtime)*365.25/tau)) - amp*np.log(1+((timespan[0]-eqtime)*365.25/tau))
    return post





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
            popt.append([-5000, 5000])
        if flag[i] == 'VELOCITY':
            if cor.correct_velo == True:
                if len(np.where(cor.velsite == site)[0]) > 0:
                    if component == 'E':
                        ve = cor.veldata[cor.velsite==site][0,2]
                        popt.append([ve-1.0, ve+2.0])
                    elif component == 'N':
                        vn = cor.veldata[cor.velsite==site][0,3]
                        popt.append([vn-1.0, vn+2.0])
                    else:
                        popt.append([-5, 5])
                else:
                    popt.append([-100, 100])
            else:
                popt.append([-100, 100])
        if flag[i] == 'EQOFFSET':
            popt.append([-5000, 5000])
        if flag[i] == 'EQDECAY':
            popt.append([-5000, 5000])
        if flag[i] == 'TAU':
            popt.append([5, 60])
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
    eqfile     = './eq_rename.cors'
    velfile    = './model.vel.gmtvec'
    offsetfile = ''
    periodfile = ''
    eq         = eqcatalog(eqfile)
    bk         = breakcatalog(eqfile)
    eqp        = eqPostList(eqfile, eq)
    cor        = correction(velfile, offsetfile, periodfile)
    mod_dict   = {'detrend':False, 'debreak': False, 'deeqoffset': False, 'depost': True, 'deseason':False}
    plot_dict  = mod_dict
    plot_dict.update({'figformat':'png', 'showfig':False})
    timespan   = [2011, 2025]
    nburns     = args.nburns
    nsteps     = args.nsteps
    post_tspan = [2021.389, 2021.631]
    eqtime     = 2021.389,
    sitefile   = 'site.list'
    sitelist   = np.genfromtxt(sitefile, dtype='S')
    
    param_dict = {
                  'constant'  : True,
                  'linear'    : True,
                  'eqlist'    : eq.eqlist,
                  'eqpostlist': eqp.eqpostlist,
                  'brklist'   : bk.breaklist,
                  'ANN'       : False,
                  'SANN'      : False,
                  'correct'   : cor}


    fid = open('postseismic.gmtvec', 'a')
    for i in range(len(sitelist)):
        posfiles = glob.glob('./'+sitelist[i].decode()+'*.pos')
        print(posfiles)
        for posfile in posfiles:
            data = posData(posfile)

            # East component
            param_dict['eqlist'] = eq.eqlist
            erun  = tsfitting(data.site, data.lon, data.lat, data.decyr, data.E, data.SE, param_dict, 'E', timespan)
            flag  = erun.flag2
            if 'EQDECAY' not in flag:
                continue
            ndim     = len(flag)
            nwalkers = 2*ndim
            popt     = set_bound(flag, cor, data.site, component = 'E')
            print(flag)
    
            starting_guess = np.random.random((nwalkers, ndim))
            for i in range(ndim):
                starting_guess[:,i] = np.random.uniform(min(popt[i]), max(popt[i]), nwalkers)
    
            sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=[popt, erun])
            sampler.run_mcmc(starting_guess, nsteps, progress=True)
            chain   = sampler.get_chain()
            np.savez('chain', chain)
            e_trace = chain[nburns:,:,:].reshape(-1, ndim)
            fig     = corner.corner(e_trace, show_titles=True)
            fig.savefig("{}_{}_posterior.png".format(data.site, "E"))
            plt.cla()
            plt.close('all')

            post_e     = postsesimic_trace(flag, e_trace, eqtime, post_tspan, method='LOG')
            cov_e      = np.mean(abs(np.diff(corner.quantile(post_e, [0.16, 0.5, 0.84]))))**2
            post_e     = corner.quantile(post_e, [0.5])[0]
            eparam     = np.array([corner.quantile(e_trace[:,i], [0.5])[0] for i in range(ndim)])
            erun.param = eparam
            erun.ifun  = erun.full_filter(erun.t) 
            erun.wrms  = np.sqrt(cov_e)
            erun.res   = erun.obs - erun.ifun(nrun.t, *eparm)


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
    
            sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=[popt, nrun])
            sampler.run_mcmc(starting_guess, nsteps, progress=True)
            chain   = sampler.get_chain()
            np.savez('{}_{}_chain'.format(data.site, "N"), chain)
            n_trace = chain[nburns:,:,:].reshape(-1, ndim)
            fig     = corner.corner(n_trace, show_titles=True)
            fig.savefig("{}_{}_posterior.png".format(data.site, "N"))

            post_n  = postsesimic_trace(flag, n_trace, eqtime, post_tspan, method='LOG')
            cov_n   = np.mean(abs(np.diff(corner.quantile(post_n, [0.16, 0.5, 0.84]))))**2
            post_n  = corner.quantile(post_n, [0.5])[0]
            nparam     = np.array([corner.quantile(n_trace[:,i], [0.5])[0] for i in range(ndim)])
            nrun.param = nparam
            nrun.ifun  = nrun.full_filter(nrun.t) 
            nrun.wrms  = np.sqrt(cov_n)
            nrun.res   = nrun.obs - nrun.ifun(nrun.t, *eparm)


            # Up component
            param_dict['eqlist'] = eq.eqlist
            urun  = tsfitting(data.site, data.lon, data.lat, data.decyr, data.U, data.SU, param_dict, 'U', timespan)
            flag  = urun.flag2
            ndim  = len(flag)
            nwalkers = 2*ndim
            popt  = set_bound(flag, cor, data.site, component = 'U')
            print(flag)
    
            starting_guess = np.random.random((nwalkers, ndim))
            for i in range(ndim):
                starting_guess[:,i] = np.random.uniform(min(popt[i]), max(popt[i]), nwalkers)
    
            sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=[popt, urun])
            sampler.run_mcmc(starting_guess, nsteps, progress=True)
            chain   = sampler.get_chain()
            np.savez('{}_{}_chain'.format(data.site, "U"), chain)
            u_trace = chain[nburns:,:,:].reshape(-1, ndim)
            fig         = corner.corner(u_trace, show_titles=True)
            fig.savefig("{}_{}_posterior.png".format(data.site, "U"))


            post_u  = postsesimic_trace(flag, u_trace, eqtime, post_tspan, method='LOG')
            cov_u   = np.mean(abs(np.diff(corner.quantile(post_u, [0.16, 0.5, 0.84]))))**2
            post_u  = corner.quantile(post_u, [0.5])[0]
            uparam     = np.array([corner.quantile(u_trace[:,i], [0.5])[0] for i in range(ndim)])
            urun.param = uparam
            urun.ifun  = urun.full_filter(urun.t) 
            urun.wrms  = np.sqrt(cov_u)
            urun.res   = urun.obs - urun.ifun(urun.t, *uparm)

            plot_obs_mod(nrun, erun, urun, nparam, eparam, uparam, plot_dict)
#           output_postseismic_ts(nrun, erun, urun, post_tspan, otype='obs')
            fid.write('{:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:5s} {:10.3f} {:10.3f}\n'.format(
                nrun.lon, nrun.lat, post_e, post_n, 
                np.sqrt(cov_e), np.sqrt(cov_n), 0.0, nrun.site, post_u, np.sqrt(cov_u)))
    fid.close()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Fit time series using MCMC method")
    parser.add_argument('--nsteps', type=str, required=False, default=10000)
    parser.add_argument('--nburns',  type=str, required=False, default=6000)
    parser.add_argument('--nwalkers', type=str, required=False, default=10)
    args = parser.parse_args()
    main(args)
