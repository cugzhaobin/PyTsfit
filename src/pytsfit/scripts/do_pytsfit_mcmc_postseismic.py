#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 19:48:51 2020

@author: zhao

Postseismic fitting with MCMC sampling. Configuration is read from a YAML
file (see config.yaml); the ``dict_fit`` section is forwarded to
``tsfitting`` as ``fit_opts`` so the quality-control options of the
refactored package apply.

For every site the postseismic displacement over ``dict_output.eqpost_tspan``
is written to ``postseismic.gmtvec``. ``--eqtime`` overrides the earthquake
epoch; by default the epoch of the first postseismic event is used.
"""

from pytsfit.data import posData, neuData
from pytsfit.models import eqcatalog, breakcatalog, eqPostList, correction
from pytsfit.tsfitting import tsfitting
from pytsfit.output import plot_obs_mod
import glob, emcee, corner, argparse, logging, yaml, os
import numpy as np
import matplotlib.pyplot as plt


def postsesimic_trace(flag, trace, eqtime, timespan, method='LOG'):
    '''
    Postseismic displacement (mm) over ``timespan`` from an MCMC trace.

    Input:
        flag     = parameter flags (flag2 of the tsfitting run)
        trace    = (nsamples, ndim) chain of posterior samples
        eqtime   = earthquake epoch in decimal year
        timespan = [t0, t1] postseismic integration window
        method   = 'LOG' or 'EXP' decay model
    Output:
        post = postseismic displacement for each sample
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
    # Points removed by the pre-fit sigma screen (fit_opts['max_sigma']) are
    # excluded from the likelihood, keeping MCMC consistent with doFitting().
    good = run.good if hasattr(run, 'good') else np.ones_like(run.obs, dtype=bool)
    res  = (obs[good]-mod[good])/run.sigma[good]
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
    logging.basicConfig(level=logging.INFO,
        format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
        datefmt="%d-%m-%Y %H:%M:%S")
    with open(args.cfgfile, 'r') as fid:
        lines = fid.read()
        cfg   = yaml.load(lines, Loader=yaml.FullLoader)
    nburns     = args.nburns
    nsteps     = args.nsteps
    dict_input = cfg['dict_input']
    eqfile     = dict_input['eqfile']
    prior_velfile    = dict_input['prior_velfile']
    prior_offsetfile = dict_input['prior_offsetfile']
    prior_periodfile = dict_input['prior_periodfile']
    sitefile   = dict_input['sitefile']
    tsdir      = dict_input['tsdir']
    tsformat   = dict_input['tsformat']
    constraint = correction(prior_velfile, prior_offsetfile, prior_periodfile)
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
    # Older config files may not have a dict_fit section; fall back to defaults.
    fit_opts   = cfg.get('dict_fit', {})
    param_dict = {
              'constant'  : dict_param['constant'],
              'linear'    : dict_param['linear'],
              'ANN'       : dict_param['annual'],
              'SANN'      : dict_param['semiannual'],
              'eqlist'    : [],
              'eqpostlist': [],
              'brklist'   : [],
              'correct'   : constraint}
    eq = None
    if dict_param['eqoffset_ne'] == True or dict_param['eqoffset_up'] == True:
        eq = eqcatalog(eqfile)
        param_dict['eqlist'] = eq.eqlist
    if dict_param['break'] == True:
        bk = breakcatalog(eqfile)
        param_dict['brklist'] = bk.breaklist
    if dict_param['eqpost_ne'] == True or dict_param['eqpost_up'] == True:
        if eq is None:
            eq = eqcatalog(eqfile)
        eqp = eqPostList(eqfile, eq)
        param_dict['eqpostlist'] = eqp.eqpostlist

    # Postseismic integration window and (optionally) the event epoch.
    post_tspan = dict_output.get('eqpost_tspan', timespan)
    if len(post_tspan) != 2:
        post_tspan = timespan

    fid = open('postseismic.gmtvec', 'a')
    for i in range(len(sitelist)):
        posfiles = glob.glob('{}/{}*.{}'.format(tsdir, sitelist[i], tsformat))
        for posfile in posfiles:
            logging.info('fitting postseismic for {}'.format(posfile))
            if tsformat == 'pos':
                data = posData(posfile)
            else:
                data = neuData(posfile)

            #
            # East component
            #
            param_dict['eqlist'] = eq.eqlist
            erun  = tsfitting(data.site, data.lon, data.lat, data.decyr, data.E, data.SE,
                              param_dict, 'E', timespan, fit_opts=fit_opts)
            flag  = erun.flag2
            if 'EQDECAY' not in flag:
                continue
            eqtime = args.eqtime if args.eqtime is not None else erun.ieqpostlist[0].eq.decyr
            ndim     = len(flag)
            nwalkers = 2*ndim
            popt     = set_bound(flag, constraint, data.site, component = 'E')
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
            erun.res   = erun.obs - erun.ifun(erun.t, *eparam)
            erun.wrms  = np.sqrt(cov_e)

            #
            # North component
            #
            param_dict['eqlist'] = eq.eqlist
            nrun  = tsfitting(data.site, data.lon, data.lat, data.decyr, data.N, data.SN,
                              param_dict, 'N', timespan, fit_opts=fit_opts)
            flag  = nrun.flag2
            ndim  = len(flag)
            nwalkers = 2*ndim
            popt  = set_bound(flag, constraint, data.site, component = 'N')
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
            nrun.res   = nrun.obs - nrun.ifun(nrun.t, *nparam)
            nrun.wrms  = np.sqrt(cov_n)

            #
            # Up component
            #
            param_dict['eqlist'] = eq.eqlist
            urun  = tsfitting(data.site, data.lon, data.lat, data.decyr, data.U, data.SU,
                              param_dict, 'U', timespan, fit_opts=fit_opts)
            flag  = urun.flag2
            ndim  = len(flag)
            nwalkers = 2*ndim
            popt  = set_bound(flag, constraint, data.site, component = 'U')
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
            urun.res   = urun.obs - urun.ifun(urun.t, *uparam)
            urun.wrms  = np.sqrt(cov_u)

            plot_obs_mod(nrun, erun, urun, nparam, eparam, uparam, dict_plot)
            fid.write('{:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:5s} {:10.3f} {:10.3f}\n'.format(
                nrun.lon, nrun.lat, post_e, post_n,
                np.sqrt(cov_e), np.sqrt(cov_n), 0.0, nrun.site, post_u, np.sqrt(cov_u)))
    fid.close()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Fit postseismic time series using MCMC method")
    parser.add_argument('--cfgfile', type=str, required=True, help='configure file in YAML format')
    parser.add_argument('--sitelist', type=str, required=False, help='This will overwrite the sitefile in the configure file.', nargs='+')
    parser.add_argument('--eqtime', type=float, required=False, default=None,
                        help='earthquake epoch (decimal year); defaults to the first postseismic event')
    parser.add_argument('--nsteps', type=int, required=False, default=10000)
    parser.add_argument('--nburns',  type=int, required=False, default=6000)
    args = parser.parse_args()
    main(args)
