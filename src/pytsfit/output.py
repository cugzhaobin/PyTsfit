#!/usr/bin/env python
# ----------------------------------------------------------
# Output and plotting functions (output_* / plot_obs_mod).
#
# Extracted from PyTsfit.py during the module split. No behaviour change.
# ----------------------------------------------------------
import os, logging
import numpy as np
import matplotlib.pyplot as plt

def output_velo(nrun, erun, urun, fid=None, fmt='GMT'):
    '''
    Ouput secular velocity

    Input:
        nrun   = an instance of class tsfitting for north component
        erun   = an instance of class tsfitting for east component
        urun   = an instance of class tsfitting for up component
        fid    = fid
        fmt    = GMT/IOS3D/GLOBK/DETAIL
    '''

    if 'VELOCITY' not in nrun.flag: return
    if 'VELOCITY' not in erun.flag: return
    if 'VELOCITY' not in urun.flag: return

    print('output velo result')
    if fid is not None:
        if fmt == 'GMT':
            fid.write(" %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %s\n"
                    %(nrun.lon, nrun.lat, erun.param[1], nrun.param[1],
                      np.sqrt(np.diag(erun.cov))[1], np.sqrt(np.diag(nrun.cov))[1], 0.0, nrun.site))
        elif fmt == 'IOS3D':
            fid.write(" %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %s %10.3f %10.3f\n"
                    %(nrun.lon, nrun.lat, erun.param[1], nrun.param[1],
                      np.sqrt(np.diag(erun.cov))[1], np.sqrt(np.diag(nrun.cov))[1], 0.0, nrun.site,
                      urun.param[1], np.sqrt(np.diag(urun.cov))[1]))
        elif fmt == 'GLOBK':
            if np.sqrt(np.diag(erun.cov))[1]<5 and np.sqrt(np.diag(nrun.cov))[1]<5:
                fid.write(" %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %8s\n"
                        %(nrun.lon, nrun.lat, erun.param[1], nrun.param[1], erun.param[1], nrun.param[1],
                         np.sqrt(np.diag(erun.cov))[1], np.sqrt(np.diag(nrun.cov))[1], 0.0,
                         urun.param[1], urun.param[1], np.sqrt(np.diag(urun.cov))[1], nrun.site))
            else:
                fid.write("#%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %8s\n"
                        %(nrun.lon, nrun.lat, erun.param[1], nrun.param[1], erun.param[1], nrun.param[1],
                          np.sqrt(np.diag(erun.cov))[1], np.sqrt(np.diag(nrun.cov))[1], 0.0,
                          urun.param[1], urun.param[1], np.sqrt(np.diag(urun.cov))[1], nrun.site))
        elif fmt == 'DETAIL':
            if np.sqrt(np.diag(erun.cov))[1]<5 and np.sqrt(np.diag(nrun.cov))[1]<5:
                fid.write(" %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %8s #%10.3f %10.3f %10.3f %d\n"
                        %(nrun.lon, nrun.lat, erun.param[1], nrun.param[1], erun.param[1], nrun.param[1],
                         np.sqrt(np.diag(erun.cov))[1], np.sqrt(np.diag(nrun.cov))[1], 0.0,
                         urun.param[1], urun.param[1], np.sqrt(np.diag(urun.cov))[1], nrun.site, min(nrun.t), max(nrun.t), max(nrun.t)-min(nrun.t), len(nrun.obs)))
            else:
                fid.write("#%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %8s #%10.3f %10.3f %10.3f %d\n"
                        %(nrun.lon, nrun.lat, erun.param[1], nrun.param[1], erun.param[1], nrun.param[1],
                          np.sqrt(np.diag(erun.cov))[1], np.sqrt(np.diag(nrun.cov))[1], 0.0,
                          urun.param[1], urun.param[1], np.sqrt(np.diag(urun.cov))[1], nrun.site, min(nrun.t), max(nrun.t), max(nrun.t)-min(nrun.t), len(nrun.obs)))

def output_postseismic_velo(nrun, erun, urun, time_span, fid=None):
    '''
    Oput postseismic velocity during time_span

    Mod by Zhao Bin, Jan 31, 2019. Fix bug when generate modeled time series using get_mod() method.
    Mod by Zhao Bin, Mar 15, 2019. discard input parameter mod_dict

    Input:
        nrun   = an instance of class tsfitting for north component
        erun   = an instance of class tsfitting for east component
        urun   = an instance of class tsfitting for up component
        plot_dict = {'detrend':True, 'debreak': True, 'deeqoffset': True, 'depost': True, 'deseason':True}
    '''

    # if postseismic term is not estimated, return
    if 'EQDECAY' not in nrun.flag: return

    # define mod_dict
    mod_dict   = {'detrend':False, 'debreak': False, 'deeqoffset': False, 'depost': True, 'deseason':False}

    # north component
    nt, nm = nrun.get_mod(time_span)
    idx    = np.where(np.logical_and(nt>time_span[0], nt<time_span[1]))[0]
    if len(idx) > 0:
        nobs_correct, nmod_correct = nrun.get_correct(mod_dict)
        n_p =  np.polyfit(nt[idx], nmod_correct[idx], 1)

    # east component
    et, nm = erun.get_mod(time_span)
    idx    = np.where(np.logical_and(et>time_span[0], et<time_span[1]))[0]
    if len(idx) > 0:
        eobs_correct, emod_correct = erun.get_correct(mod_dict)
        e_p =  np.polyfit(et[idx], emod_correct[idx], 1)

    # print the results
    if len(idx) > 0:
        if fid is None:
            print("%s %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %s\n"
                %(" VEL ", nrun.lon, nrun.lat, e_p[0], n_p[0], 0.0, 0.0, 0.0, nrun.site))
        else:
            fid.write(" %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %s\n"
                    %(nrun.lon, nrun.lat, e_p[0], n_p[0], erun.wrms/2.0, nrun.wrms/2.0, 0.0, nrun.site))


def output_postseismic_disp(nrun, erun, urun, time_span, fid=None, eqcode=None):
    '''
    Oput postseismic displacement during time_span

    Mod by Zhao Bin, Jan. 31, 2019. Fix bug when generate modeled time series using get_mod() method.
    Mod by Zhao Bin. Feb. 17, 2019. Ouput vertical displacements.
    Mod by Zhao Bin, Mar. 15, 2019. discard input parameter mod_dict

    Input:
        nrun   = an instance of class tsfitting for north component
        erun   = an instance of class tsfitting for east component
        urun   = an instance of class tsfitting for up component
        time_span = time range for postseismic calculation
        fid    = file ID
        eqcode = earthquake code, only postseismic displacement associated with this EQ will be output
    '''

    # if postseismic term is not estimated, return
    if 'EQDECAY' not in nrun.flag: return

    # define mod_dict
    mod_dict   = {'detrend':False, 'debreak': False, 'deeqoffset': False, 'depost': True, 'deseason':False}

    # north component
    nt, nm = nrun.get_mod()
    idx    = np.where(np.logical_and(nt>time_span[0], nt<time_span[1]))[0]
    if len(idx) > 0:
        nobs_correct, nmod_correct = nrun.get_correct(mod_dict, time_span=time_span, eqcode=eqcode)
        if len(nmod_correct) == 0:
            logging.warning("No correction model for {0:s}".format(nrun.site))
            return
        ndisp = nmod_correct[-1]-nmod_correct[0]

    # east component
    et, em = erun.get_mod()
    idx    = np.where(np.logical_and(et>time_span[0], et<time_span[1]))[0]
    if len(idx) > 0:
        eobs_correct, emod_correct = erun.get_correct(mod_dict, time_span=time_span, eqcode=eqcode)
        if len(emod_correct) == 0:
            logging.warning("No correction model for {0:s}".format(erun.site))
            return
        edisp = emod_correct[-1]-emod_correct[0]

    # vertical component
    ut, um = urun.get_mod()
    idx    = np.where(np.logical_and(ut>time_span[0], ut<time_span[1]))[0]
    if len(idx) > 0:
        uobs_correct, umod_correct = urun.get_correct(mod_dict, time_span=time_span, eqcode=eqcode)
        if len(umod_correct) == 0:
            logging.warning("No correction model for {0:s}".format(urun.site))
            return
        udisp = umod_correct[-1]-umod_correct[0]

    # print the results
    if len(idx) > 0:
        if fid is None:
            print("%s %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %s\n"
                %(" VEL ", nrun.lon, nrun.lat, edisp, ndisp, 0.0, 0.0, 0.0, nrun.site))
        else:
            fid.write(" %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %s %10.3f %10.3f\n" \
                    %(nrun.lon, nrun.lat, edisp, ndisp, erun.wrms, nrun.wrms, 0.0, nrun.site,\
                      udisp, urun.wrms))

def output_eqoffset(nrun, erun, urun, fid=None, fmt='GMT2D'):
    '''
    Output coseismic displacements.
    Written by Zhao Bin, Feb. 21, 2019.
    Mod by Zhao Bin, Mar. 28, 2019. Fix bug when no EQOFFSET in flag2
    The error is very large

    Input:
        nrun   = an instance of class tsfitting for north component
        erun   = an instance of class tsfitting for east  component
        urun   = an instance of class tsfitting for up    component
        fid    = file ID to output
        fmt    = format of output file. GMT2D/GMT3D/IOS3D

    '''
    if 'EQOFFSET' not in nrun.flag2: return
    idx = nrun.flag2.tolist().index('EQOFFSET')
    for i in range(len(nrun.ieqlist)):
        if fid is not None:
            if fmt == 'GMT2D':
                fid.write('# earthquake code = %2s, decyr = %10.3f\n' %(nrun.ieqlist[i].code, nrun.ieqlist[i].decyr))
                fid.write("%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %s # %s\n"
                    %(nrun.lon, nrun.lat, erun.param[idx+i], nrun.param[idx+i], np.sqrt(np.diag(erun.cov))[idx+i],
                    np.sqrt(np.diag(nrun.cov))[idx+i], 0.0, nrun.site, nrun.ieqlist[i].code))
            if fmt == 'GMT3D':
                fid.write('# earthquake code = %2s, decyr = %10.3f\n' %(nrun.ieqlist[i].code, nrun.ieqlist[i].decyr))
                if 'EQOFFSET' not in urun.flag2:
                    logging.info(' No offsets due to any earthquakes are estimated.')
                    continue
                fid.write("%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %s # %s\n"
                    %(nrun.lon, nrun.lat, erun.param[idx+i], nrun.param[idx+i], urun.param[idx+i], np.sqrt(np.diag(erun.cov))[idx+i],
                    np.sqrt(np.diag(nrun.cov))[idx+i], np.sqrt(np.diag(urun.cov))[idx+i], nrun.site, nrun.ieqlist[i].code))
            if fmt == 'IOS3D':
                fid.write('# earthquake code = %2s, decyr = %10.3f\n' %(nrun.ieqlist[i].code, nrun.ieqlist[i].decyr))
                if 'EQOFFSET' not in urun.flag2:
                    logging.info('No offsets due to any earthquakes are estimated.')
                    continue
                fid.write("%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %s %10.3f %10.3f # %s\n"
                    %(nrun.lon, nrun.lat, erun.param[idx+i], nrun.param[idx+i], np.sqrt(np.diag(erun.cov))[idx+i],
                    np.sqrt(np.diag(nrun.cov))[idx+i], 0.0, nrun.site, urun.param[idx+i], np.sqrt(np.diag(urun.cov))[idx+i],
                    nrun.ieqlist[i].code))

def output_break(nrun, erun, urun, fid=None):
    '''
    Ouput non-earthquake breaks.

    Written by Zhao Bin, Feb. 21, 2019.

    Input:
        nrun   = an instance of class tsfitting for north component
        erun   = an instance of class tsfitting for east  component
        urun   = an instance of class tsfitting for up    component
    '''
    idx = nrun.flag2.tolist().index('BREAK')
    if fid is not None:
        for i in range(len(nrun.ibrklist)):
            print("%10.3f %10.3f %10.3f" %(erun.param[idx+i], nrun.param[idx+i], nrun.ibrklist[i].decyr))

def output_postseismic_ts(nrun, erun, urun, time_span, otype='obs', eqcode=None):
    '''
    Output postseismic time series during time_span.

    Written by Zhao Bin, Mar. 15, 2019.

    Input:
        nrun      = an instance of class tsfitting for north component
        erun      = an instance of class tsfitting for east  component
        urun      = an instance of class tsfitting for up    component
        time_span = time range for output
        otype     = obs/mod to ouput observed/modeled time series (continuous)
        eqcode    = earthquake code
    '''
    # if postseismic term is not estimated, return
    if 'EQDECAY' not in nrun.flag: return

    # define mod_dict
    mod_dict   = {'detrend': True, 'debreak': True, 'deeqoffset': True, 'depost': False, 'deseason': True}
    nobs_correct, _ = nrun.get_correct(mod_dict, eqcode=eqcode)
    eobs_correct, _ = erun.get_correct(mod_dict, eqcode=eqcode)
    uobs_correct, _ = urun.get_correct(mod_dict, eqcode=eqcode)

    # pure postseismic
    idx  = np.where(np.logical_and(nrun.t>time_span[0], nrun.t<time_span[1]))[0]
    nobs = (nrun.obs-nobs_correct)/1e3
    eobs = (erun.obs-eobs_correct)/1e3
    uobs = (urun.obs-uobs_correct)/1e3

    # output to file
    if eqcode != None:
        outfile = f'{nrun.site}_{eqcode}.neu'
    else:
        outfile = f'{nrun.site}_pos.neu'
    with open(outfile, 'w') as fid:
        fid.write('# YEAR_SINCE_EQ  NORTH(m)  EAST(m) Up(m)  Nsigma(m)  Esigma(m)  Usigma(m)\n')
        if otype == 'mod':
            pass
        elif otype == 'obs':
            n0, e0, u0 = nobs[idx[0]], eobs[idx[0]], uobs[idx[0]]
            for i in idx:
                fid.write("%12.5f %10.4f %10.4f %10.4f %10.4f  %10.4f %10.4f\n"
                        %(nrun.t[i]-time_span[0], nobs[i]-n0, eobs[i]-e0, uobs[i]-u0,
                            nrun.sigma[i]/1e3, erun.sigma[i]/1e3, urun.sigma[i]/1e3))


def output_period(nrun, erun, urun, nparam, eparam, uparam):
    '''
    Output estimation of period term

    Written by Zhao Bin, Mar.  9, 2020.

    Input:
        nrun   = an instance of class tsfitting for north component
        erun   = an instance of class tsfitting for east  component
        urun   = an instance of class tsfitting for up    component
        nparam = estimated parameters for north component
        eparam = estimated parameters for east  component
        uparam = estimated parameters for up    component
    '''

    period = np.zeros(12)

    # east component
    idx    =  np.where(erun.flag2 == 'ANNUAL_SIN')[0]
    if len(idx)>0:
        period[0] = eparam[idx]
        period[1] = eparam[idx+1]
    idx    =  np.where(erun.flag2 == 'SANNUAL_SIN')[0]
    if len(idx)>0:
        period[2] = eparam[idx]
        period[3] = eparam[idx+1]

    # north component
    idx    =  np.where(nrun.flag2 == 'ANNUAL_SIN')[0]
    if len(idx)>0:
        period[4] = nparam[idx]
        period[5] = nparam[idx+1]
    idx    =  np.where(nrun.flag2 == 'SANNUAL_SIN')[0]
    if len(idx)>0:
        period[6] = nparam[idx]
        period[7] = nparam[idx+1]

    # up component
    idx    =  np.where(urun.flag2 == 'ANNUAL_SIN')[0]
    if len(idx)>0:
        period[8] = uparam[idx]
        period[9] = uparam[idx+1]
    idx    =  np.where(urun.flag2 == 'SANNUAL_SIN')[0]
    if len(idx)>0:
        period[10] = uparam[idx]
        period[11] = uparam[idx+1]

    with open('period.dat', 'a') as fid:
        fid.write("{:4.1f} {:4.1f} {:4.1f} {:4.1f} {:4.1f} {:4.1f} {:4.1f} {:4.1f} {:4.1f} {:4.1f} {:4.1f} {:4.1f} {:s}\n".format(
            period[0], period[1], period[2], period[3], period[4], period[5], period[6],
            period[7], period[8], period[9], period[10], period[11], nrun.site))

def output_summary(nrun, erun, urun, fid=None):
    '''
    Oput summary of tsfitting

    Write by Zhao Bin, Oct. 11, 2022.

    Input:
        nrun   = an instance of class tsfitting for north component
        erun   = an instance of class tsfitting for east  component
        urun   = an instance of class tsfitting for up    component
    '''

    fid.write(" %10.3f %10.3f %10.3f %10.3f %10.3f %10d %10.3f %10.3f %10.3f %6s\n"
                %(nrun.lon, nrun.lat, min(nrun.t), max(nrun.t), max(nrun.t)-min(nrun.t), len(nrun.obs), erun.wrms, nrun.wrms, urun.wrms, nrun.site))



def plot_obs_mod(nrun, erun, urun, nparam, eparam, uparam, plot_dict, nwrms=5, nsigma=3):
    '''
    Plot observed and modeled time series.
    Mod by Zhao Bin, Jan. 10, 2019. Fix bug when to estimation is done.
    Mod by Zhao Bin, Jan. 31, 2019. Plot error bar and labels
    Mod by Zhao Bin, Feb. 17, 2019. Ignore large uncertainities
    Mod by Zhao Bin, Apr. 27, 2023. plot figure in html format

    Input:
        nrun   = an instance of class tsfitting for north component
        erun   = an instance of class tsfitting for east  component
        urun   = an instance of class tsfitting for up    component
        nparam = estimated parameters for north component
        eparam = estimated parameters for east  component
        uparam = estimated parameters for up    component
        plot_dict = {'detrend':True, 'debreak': True, 'deeqoffset': True, 'depost': True, 'deseason':True,
                     'figformat': 'jpg', 'showfig': True}
    '''
    if len(nrun.flag) == 0: return

    #####################################################################
    # North component
    #####################################################################
    nt, nm = nrun.get_mod()
    nobs_correct, nmod_correct = nrun.get_correct(plot_dict)
    pidx    = np.where(np.logical_and(nrun.sigma<nwrms*nrun.wrms, abs(nrun.res)<nsigma*nrun.sigma))[0]
    north_t = nrun.t[pidx]
    north_d = nrun.obs[pidx]-nobs_correct[pidx]
    north_e = nrun.sigma[pidx]
    nidx    =  np.where(nrun.flag2 == 'VELOCITY')[0]


    #####################################################################
    # East component
    #####################################################################
    et, em = erun.get_mod()
    eobs_correct, emod_correct = erun.get_correct(plot_dict)
    pidx   = np.where(np.logical_and(erun.sigma<nwrms*erun.wrms, abs(erun.res)<nsigma*erun.sigma))[0]
    east_t = erun.t[pidx]
    east_d = erun.obs[pidx]-eobs_correct[pidx]
    east_e = erun.sigma[pidx]
    eidx   = np.where(erun.flag2 == 'VELOCITY')[0]


    #####################################################################
    # Vertical component
    #####################################################################
    ut, um = urun.get_mod()
    uobs_correct, umod_correct = urun.get_correct(plot_dict)
    pidx = np.where(np.logical_and(urun.sigma<nwrms*urun.wrms, abs(urun.res)<nsigma*urun.sigma))[0]
    up_t = urun.t[pidx]
    up_d = urun.obs[pidx]-uobs_correct[pidx]
    up_e = urun.sigma[pidx]
    uidx = np.where(urun.flag2 == 'VELOCITY')[0]

    if plot_dict['figformat'] == 'html':
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        fig = make_subplots(rows=3, cols=1,
            subplot_titles=('North position timeseries', 'East position timeseries', 'Vertical position timeseries'))
        fig.add_trace(go.Scatter(x=north_t, y=north_d, mode='markers', marker={'size': 13}, name='N'), row=1, col=1)
        fig.add_trace(go.Scatter(x=east_t, y=east_d, mode='markers', marker={'size': 13}, name='E'), row=2, col=1)
        fig.add_trace(go.Scatter(x=up_t, y=up_d, mode='markers', marker={'size': 13}, name='U'), row=3, col=1)
        fig.add_trace(go.Scatter(x=nt, y=nm-nmod_correct, mode='lines', line={'width': 3}, name='N'), row=1, col=1)
        fig.add_trace(go.Scatter(x=et, y=em-emod_correct, mode='lines', line={'width': 3}, name='E'), row=2, col=1)
        fig.add_trace(go.Scatter(x=ut, y=um-umod_correct, mode='lines', line={'width': 3}, name='U'), row=3, col=1)
        fig.update_xaxes(title_text='Year',row=1, col=1)
        fig.update_xaxes(title_text='Year',row=2, col=1)
        fig.update_xaxes(title_text='Year',row=3, col=1)
        fig.update_yaxes(title_text='North',row=1, col=1)
        fig.update_yaxes(title_text='East',row=2, col=1)
        fig.update_yaxes(title_text='Up',row=3, col=1)
        fig.update_layout(autosize=True, height=1500, title='Position time series at {}'.format(nrun.site))
        with open('{}.html'.format(nrun.site), 'w') as f:
            f.write(fig.to_html(full_html=False, include_plotlyjs='directory'))
        if plot_dict['showfig']:
            fig.show()
    else:
        plt.figure(figsize=(9,12))
        plt.subplots_adjust(top=0.8)
        plt.subplot(3,1,1)
        if len(nidx) == 0:
            plt.title(r'vel=%.2f $\pm$ %.2f mm/yr' %(0.0, 0.0), loc='left')
        else:
            plt.title(r'vel=%.2f $\pm$ %.2f mm/yr' %(nparam[1], 0.0), loc='left')
#           plt.title(r'vel=%.2f $\pm$ %.2f mm/yr' %(nparam[1], np.sqrt(np.diag(nrun.cov))[1]), loc='left')

        plt.errorbar(north_t, north_d, yerr=north_e, ecolor='gray',
            elinewidth=0.2, capsize=1, capthick=0.5, fmt='o', ms=3, mfc='r', mec='black', mew=0, zorder=1)
        plt.plot(nt, nm-nmod_correct, color='black', zorder=2)
        plt.ylabel('North (mm)')

        plt.subplot(3,1,2)
        if len(eidx) == 0:
            plt.title(r'vel=%.2f $\pm$ %.2f mm/yr' %(0.0, 0.0), loc='left')
        else:
            plt.title(r'vel=%.2f $\pm$ %.2f mm/yr' %(eparam[1], 0.0), loc='left')
#           plt.title(r'vel=%.2f $\pm$ %.2f mm/yr' %(eparam[1], np.sqrt(np.diag(erun.cov))[1]), loc='left')

        plt.errorbar(east_t, east_d, yerr=east_e, ecolor='gray',
            elinewidth=0.2, capsize=1, capthick=0.5, fmt='o', ms=3, mfc='g', mec='black', mew=0, zorder=1)
        plt.plot(et, em-emod_correct, color='black', zorder=2)
        plt.ylabel('East (mm)')


        plt.subplot(3,1,3)
        if len(uidx) == 0:
            plt.title(r'vel=%.2f $\pm$ %.2f mm/yr' %(0.0, 0.0), loc='left')
        else:
            plt.title(r'vel=%.2f $\pm$ %.2f mm/yr' %(uparam[1], 0.0), loc='left')
#           plt.title(r'vel=%.2f $\pm$ %.2f mm/yr' %(uparam[uidx], np.sqrt(np.diag(urun.cov))[1]), loc='left')
        plt.errorbar(up_t, up_d, yerr=up_e, ecolor='gray',
            elinewidth=0.2, capsize=1, capthick=0.5, fmt='o', ms=3, mfc='b', mec='black', mew=0, zorder=1)
        plt.plot(ut, um-umod_correct, color='black', zorder=2)
        plt.ylabel('Up (mm)')
        plt.xlabel('Time (year)')
        plt.suptitle(nrun.site)

        # Adjust the space
        plt.subplots_adjust(hspace=0.6, wspace=0.6)

        # plot title
        plt.suptitle("Time Series of Site Position drawn by zhao at Institute of Seismology\n\nStation: "+nrun.site+"\n\n %10.3fN %10.3fE %6.2f(m)\n\n %d Daily solution (%7.2f-%7.2f)" %(nrun.lat, nrun.lon, 0.0, len(nrun.t), min(nt), max(nt)), fontsize=15)
#   plt.tight_layout()

        # show the figure
        if plot_dict['showfig']:
            plt.show()
        fmt = plot_dict['figformat']
        plt.savefig(nrun.site+'.'+fmt, format=fmt, dpi=300)
        plt.close()

def output_obs_mod(nrun, erun, urun, nparam, eparam, uparam, mod_dict, obs=True, mod=True):
    '''
    Output observed and modeled time series.

    Input:
        nrun   = an instance of class tsfitting for north component
        erun   = an instance of class tsfitting for east  component
        urun   = an instance of class tsfitting for up    component
        nparam = estimated parameters for north component
        eparam = estimated parameters for east  component
        uparam = estimated parameters for up    component
        mod_dict = {'detrend':True, 'debreak': True, 'deeqoffset': True, 'depost': True, 'deseason':True}
    '''
    if obs==False and mod==False: return
    if len(nrun.flag) == 0: return

    #####################################################################
    # North component
    #####################################################################
    nt, nm = nrun.get_mod()
    idx    =  np.where(nrun.flag2 == 'VELOCITY')[0]
    nobs_correct, nmod_correct = nrun.get_correct(mod_dict)

    #####################################################################
    # East component
    #####################################################################
    _, em = erun.get_mod()
    idx    =  np.where(erun.flag2 == 'VELOCITY')[0]
    eobs_correct, emod_correct = erun.get_correct(mod_dict)

    #####################################################################
    # Vertical component
    #####################################################################
    _, um = urun.get_mod()
    idx    =  np.where(urun.flag2 == 'VELOCITY')[0]
    uobs_correct, umod_correct = urun.get_correct(mod_dict)

    if obs==True:
        fname = '{}_obs.dat'.format(nrun.site)
        with open(fname, 'w') as fid:
            fid.write("# Reference location: {:12.3f} {:12.3f}\n".format(nrun.lon, nrun.lat))
            fid.write("#      Decyr        N(mm)        E(mm)        U(mm)       Sn(mm)       Se(mm)       Su(mm)   Nflag   Eflag   Uflag\n")

        with open(fname, 'a') as fid:
            flag = np.zeros((len(nrun.t), 3), dtype=int)
            idx  = np.where(abs(nrun.res)>3*nrun.wrms)[0]
            flag[idx,0] = 1
            idx  = np.where(abs(erun.res)>3*erun.wrms)[0]
            flag[idx,1] = 1
            idx  = np.where(abs(urun.res)>3*urun.wrms)[0]
            flag[idx,2] = 1

            idx  = np.where(nrun.sigma>20)[0]
            flag[idx,0] = 2
            idx  = np.where(erun.sigma>20)[0]
            flag[idx,1] = 2
            idx  = np.where(urun.sigma>50)[0]
            flag[idx,2] = 2

            for i in range(len(nrun.t)):
                fid.write("{:12.5f} {:12.2f} {:12.2f} {:12.2f} {:12.2f} {:12.2f} {:12.2f} {:7d} {:7d} {:7d}\n".format(nrun.t[i],
                           nrun.obs[i]-nobs_correct[i],
                           erun.obs[i]-eobs_correct[i],
                           urun.obs[i]-uobs_correct[i],
                           nrun.sigma[i],
                           erun.sigma[i],
                           urun.sigma[i],
                           flag[i,0], flag[i,1], flag[i,2]))

    if mod==True:
        fname = '{}_mod.dat'.format(nrun.site)
        with open(fname, 'w') as fid:
            fid.write("# {} {:12.2f} {:12.2f}\n".format(nrun.site, nrun.lon, nrun.lat))
            fid.write("#      Decyr        N(mm)        E(mm)        U(mm)\n")

        with open(fname, 'a') as fid:
            for i in range(len(nt)):
                fid.write("{:12.5f} {:12.2f} {:12.2f} {:12.2f}\n".format(nt[i],
                           nm[i]-nmod_correct[i],
                           em[i]-emod_correct[i],
                           um[i]-umod_correct[i]))

def output_param(nrun, erun, urun, nparam, eparam, uparam):
    '''
    Output residual files like tsfit program.

    Input:
        nrun   = an instance of class tsfitting for north component
        erun   = an instance of class tsfitting for east  component
        urun   = an instance of class tsfitting for up    component
        nparam = estimated parameters for north component
        eparam = estimated parameters for east  component
        uparam = estimated parameters for up    component
    '''

    parm = {}
    parm['site'] = nrun.site
    idx          = np.where(nrun.flag2 == 'VELOCITY')[0]
    parm['Vn']   = nparam[idx]
    parm['Sn']   = np.sqrt(np.diag(nrun.cov))[idx]
    idx          = np.where(erun.flag2 == 'VELOCITY')[0]
    parm['Ve']   = eparam[idx]
    parm['Se']   = np.sqrt(np.diag(erun.cov))[idx]
    idx          = np.where(urun.flag2 == 'VELOCITY')[0]
    parm['Vu']   = uparam[idx]
    parm['Su']   = np.sqrt(np.diag(urun.cov))[idx]
    if 'EQOFFSET' not in nrun.flag2:
        parm['eqoffset'] = None
        parm['eqinfo']   = None
    else:
        idx = nrun.flag2.tolist().index('EQOFFSET')
        eqcode = []
        eqyear = []
        eqdisp = []
        for i in range(len(nrun.ieqlist)):
            eqcode.append(nrun.ieqlist[i].code)
            eqyear.append(nrun.ieqlist[i].decyr)
            idisp = [erun.param[idx+i], nrun.param[idx+i], urun.param[idx+i], \
                    np.sqrt(np.diag(erun.cov))[idx+i], \
                    np.sqrt(np.diag(nrun.cov))[idx+i], \
                    np.sqrt(np.diag(urun.cov))[idx+i]]
            eqdisp.append(idisp)


        eqinfo = {}
        eqinfo['eqcode'] = eqcode
        eqinfo['eqyear'] = eqyear
        eqinfo['eqdisp'] = eqdisp
        parm['eqinfo']   = eqinfo
    np.save('{}_par'.format(nrun.site), parm)
