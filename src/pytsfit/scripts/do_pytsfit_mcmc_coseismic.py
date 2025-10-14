#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 19:48:51 2020

@author: zhao
"""

from PyTsfit import posData, correction, tsfitting, neuData
from PyTsfit import eqcatalog, breakcatalog, eqPostList
from PyTsfit import plot_obs_mod, output_param, output_obs_mod,output_velo,output_eqoffset
from PyTsfit import output_break,output_postseismic_disp,output_postseismic_ts
import glob, emcee, corner, argparse, logging, yaml, os
import numpy as np
import matplotlib.pyplot as plt


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
    with open(args.cfgfile, 'r') as fid:
        lines = fid.read()
        cfg   = yaml.load(lines, Loader=yaml.FullLoader)
    nburns     = args.nburns
    nsteps     = args.nsteps
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
              'correct'   : constraint}
    if args.annual == "True" or args.annual == "true":
        param_dict['ANN'] = True
    else:
        param_dict['ANN'] = False
    if args.semiannual == "True" or args.semiannual == "true":
        param_dict['SANN'] = True
    else:
        param_dict['SANN'] = False
    if dict_param['eqoffset_ne'] == True or dict_param['eqoffset_up'] == True:
        eq = eqcatalog(eqfile)
        param_dict['eqlist'] = eq.eqlist
    if dict_param['break'] == True:
        bk = breakcatalog(eqfile)
        param_dict['brklist'] = bk.breaklist
    if dict_param['eqpost_ne'] == True or dict_param['eqpost_up'] == True:
        eqp = eqPostList(eqfile, eq)
        param_dict['eqpostlist'] = eqp.eqpostlist
        
    fid_velo      = None
    fid_eq        = None
    fid_brk       = None
    fid_post_disp = None
    if len(dict_output['velfile'])>0:
        fid_velo = open(dict_output['velfile'], 'a')
    if len(dict_output['eqoffset'])>0:
        fid_eq   = open(dict_output['eqoffset'], 'a')
    if len(dict_output['break'])>0:
        fid_brk  = open(dict_output['break'], 'a')
    if len(dict_output['eqpostdisp']) > 0:
        fid_post_disp = open(dict_output['eqpostdisp'], 'a')
    
    for i in range(len(sitelist)):
        posfiles = glob.glob('{}/{}*.{}'.format(tsdir, sitelist[i], tsformat))
        for j in range(len(posfiles)):
            logging.info('fitting time series for {}'.format(posfiles[j]))
            posfile = posfiles[j]
            if tsformat == 'pos':
                data    = posData(posfile)
            if tsformat == 'neu':
                data    = neuData(posfile)

            #
            # East component
            #
            if dict_param['eqoffset_ne'] == False: param_dict['eqlist']=[]
            if dict_param['eqpost_ne'] == False: param_dict['eqpostlist']=[]
            erun           = tsfitting(data.site, data.lon, data.lat, data.decyr, data.E, data.SE, param_dict, 'E', timespan)
            flag           = erun.flag2
            ndim           = len(flag)
            nwalkers       = 2*ndim
            popt           = set_bound(flag, constraint, data.site, component = 'E')
            starting_guess = np.random.random((nwalkers, ndim))
            for i in range(ndim):
                starting_guess[:,i] = np.random.uniform(min(popt[i]), max(popt[i]), nwalkers)
            sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=[popt, erun])
            sampler.run_mcmc(starting_guess, nsteps, progress=True)
            chain   = sampler.get_chain()
            e_trace = chain[nburns:,:,:].reshape(-1, ndim)
            eparam  = np.array([corner.quantile(e_trace[:,i], [0.5])[0] for i in range(ndim)])
            erun.param = eparam
            erun.ifun  = erun.full_filter(erun.t) 
            cov        = np.zeros((ndim, ndim))
            for i in range(ndim):
                cov[i,i] = np.mean(abs(np.diff(corner.quantile(e_trace[:,i], [0.16, 0.5, 0.84]))))**2
            erun.cov  = cov
            erun.wrms = np.sqrt(cov)
            fig       = corner.corner(e_trace[:,[0,2]], show_titles=True)
            fig.savefig("{}_{}_posterior.png".format(data.site, "E"))
            plt.cla()
            plt.close('all')
            np.savez('chain', chain)

            #
            # North component
            #
            if dict_param['eqoffset_ne'] == False: param_dict['eqlist']=[]
            if dict_param['eqpost_ne'] == False: param_dict['eqpostlist']=[]
            nrun           = tsfitting(data.site, data.lon, data.lat, data.decyr, data.N, data.SN, param_dict, 'N', timespan)
            flag           = nrun.flag2
            ndim           = len(flag)
            nwalkers       = 2*ndim
            popt           = set_bound(flag, constraint, data.site, component = 'N')
            starting_guess = np.random.random((nwalkers, ndim))
            for i in range(ndim):
                starting_guess[:,i] = np.random.uniform(min(popt[i]), max(popt[i]), nwalkers)
            sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=[popt, nrun])
            sampler.run_mcmc(starting_guess, nsteps, progress=True)
            chain   = sampler.get_chain()
            n_trace = chain[nburns:,:,:].reshape(-1, ndim)
            nparam  = np.array([corner.quantile(n_trace[:,i], [0.5])[0] for i in range(ndim)])
            nrun.param = nparam
            nrun.ifun  = nrun.full_filter(nrun.t) 
            cov        = np.zeros((ndim, ndim))
            for i in range(ndim):
                cov[i,i] = np.mean(abs(np.diff(corner.quantile(n_trace[:,i], [0.16, 0.5, 0.84]))))**2
            nrun.cov   = cov
            nrun.wrms  = np.sqrt(cov)
            nrun.res   = nrun.obs-nrun.ifun(nrun.t, *nrun.param)
            fig        = corner.corner(n_trace, show_titles=True)
            fig.savefig("{}_{}_posterior.png".format(data.site, "N"))
            np.savez('{}_{}_chain'.format(data.site, "N"), chain)

            #
            # Up component
            #
            if dict_param['eqoffset_up'] == False: param_dict['eqlist']=[]
            if dict_param['eqpost_up'] == False: param_dict['eqpostlist']=[]
            urun  = tsfitting(data.site, data.lon, data.lat, data.decyr, data.U, data.SU, param_dict, 'U', timespan)
            flag  = urun.flag2
            ndim  = len(flag)
            nwalkers = 2*ndim
            popt  = set_bound(flag, constraint, data.site, component = 'U')    
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
            uparam = np.array([corner.quantile(u_trace[:,i], [0.5])[0] for i in range(ndim)])
            urun.param = uparam
            urun.ifun = urun.full_filter(urun.t) 
            cov = np.zeros((ndim, ndim))
            for i in range(ndim):
                cov[i,i] = np.mean(abs(np.diff(corner.quantile(u_trace[:,i], [0.16, 0.5, 0.84]))))**2
            urun.cov = cov
            urun.wrms= np.sqrt(cov)


            if dict_output['tsfig'] == True:
                plot_obs_mod(nrun, erun, urun, nparam, eparam, uparam, dict_plot, nsigma=100, nwrms=100)
            if dict_output['param'] == True:
                output_param(nrun, erun, urun, nparam, eparam, uparam)
            if dict_output['obsmod'] == True:
                output_obs_mod(nrun, erun, urun, nparam, eparam, uparam, dict_plot)
            if len(dict_output['velfile']) > 0:
                output_velo(nrun, erun, urun, fid=fid_velo, fmt='GLOBK')
            if len(dict_output['eqoffset']) > 0:
                output_eqoffset(nrun, erun, urun, fid=fid_eq)
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
    return



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Fit time series using MCMC method")
    parser.add_argument('--cfgfile', type=str, required=True, help='configure file in YAML format')
    parser.add_argument('--sitelist', type=str, required=False, help='This will overwrite the sitefile in the configure file.', nargs='+')
    parser.add_argument('--showfig', type=str, required=False, help='Show figure.')
    parser.add_argument('--annual', type=str, required=False, help='Estimate annual term.')
    parser.add_argument('--semiannual', type=str, required=False, help='Estimate semiannual term.')
    parser.add_argument('--nsteps', type=int, required=False, default=4000)
    parser.add_argument('--nburns',  type=int, required=False, default=2000)
    parser.add_argument('--nwalkers', type=int, required=False, default=10)
    args = parser.parse_args()
    main(args)
