#!/usr/bin/env python
# ----------------------------------------------------------
# Fitting GPS coordinate time series or baseline time series
#
# Zhao Bin, Institute of Seismologym CEA.
# Nov 12, 2018
# ----------------------------------------------------------
import os, sys, logging
import numpy as np
import geotools as gt
import GPSTime as gpstime
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

logging.basicConfig(level=logging.INFO,
    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
    datefmt="%d-%m-%Y %H:%M:%S")

class neuData(object):
    def __init__(self, neufile):
        dat        = np.genfromtxt(neufile)
        self.decyr = dat[:,0]
        self.N     = dat[:,1]
        self.E     = dat[:,2]
        self.U     = dat[:,3]
        self.SN    = dat[:,4]
        self.SE    = dat[:,5]
        self.SU    = dat[:,6]
        
        self.lat   = float(os.popen('grep "NEU" '+neufile+"| head -n1").read().split()[6])
        self.lon   = float(os.popen('grep "NEU" '+neufile+"| head -n1").read().split()[7])
        self.site  = os.popen('grep "NEU" '+neufile+"| head -n1").read().split()[5]
        print(self.lat, self.lon, self.site)
        logging.info('Loading neuData finished.')

class posData(object):
    '''
    POSDATA is a class representing a PBO POS file.
    '''

    
    def __init__(self, posfile):
        '''
        Constructor.
        Mod by Zhao Bin, Jan. 10, 2019. Fix bug when reading pos file
        Mod by Zhao Bin, Jan. 11, 2019. Fix bug when reading pos file, use delimiter to foramt the file
        
        Input:
            posfile    = file name of a pos file
        '''
        
        # station ID in 4 char
        self.site = ""
        
        # MJD
        self.mjd  = []
        
        # Decimal yar
        self.decyr= []
        
        # North displacement in meter
        self.N    = []
        
        # East displacement in meter
        self.E    = []
        
        # Up displacement in meter
        self.U    = []
        
        # North uncertainty in meter
        self.SN   = []
        
        # East uncertainty in meter
        self.SE   = []
        
        # Up uncertainty in meter
        self.SU   = []

        # check the input file exists
        if os.path.isfile(posfile) == False:
            logging.fatal(' The time series file %s does not exist' %(posfile))
            sys.exit()
            
        # open the file and read the header
        with open(posfile) as fid:
            for line in fid:
                if line[0] != '':
                    if line.find('ID') != -1:
                        self.site = line[16:20]
                    if line.find('NEU') == 0:
                        self.lat = float(line.split()[4])
                        self.lon = float(line.split()[5])
                        self.hei = float(line.split()[6].replace('*','0'))
                elif line[0] == ' ':
                    continue
        
        # read the time series and convert to millimeter
        data     = np.genfromtxt(posfile, skip_header=37, \
                   delimiter=(9,7,11,15,15,15,9,9,9,7,7,7,19,16,11,12,10,10,11,9,9,7,7,7,6))
        if data.ndim == 1:
            data = data.reshape((1, len(data)))
        idx      = list(set(range(len(data))) - set(np.where(np.isnan(data[:,[15,16,17]]))[0]))
        self.MJD = data[idx,2]
        
        # convert to mm
        self.N   = data[idx,15]*1e3
        self.E   = data[idx,16]*1e3
        self.U   = data[idx,17]*1e3
        self.SN  = data[idx,18]*1e3
        self.SE  = data[idx,19]*1e3
        self.SU  = data[idx,20]*1e3
        idx      = np.where(self.SN == 0.0)[0]
        self.SN[idx] = 1000.0
        self.SE[idx] = 1000.0
        self.SU[idx] = 1000.0
        
        # convert the MJD to decimal year
        self.decyr = np.array([gpstime.jd_to_decyrs(self.MJD[i]) 
                                for i in range(len(self.MJD))])


    def plot_pos(self, time_range=[], show=False):
        '''
        Plot raw POS time series.
        
        Input:
            time_range = [start_time, end_time] in decimal year
        '''

        if len(time_range) == 2:
            idx = np.where(np.logical_and(self.decyr>time_range[0], 
                                          self.decyr<time_range[1]))[0]
        else:
            idx = np.arange(0,len(self.decyr))
                                          
                                      
        # North component
        plt.figure(figsize=(9,12))
        plt.subplot(3,1,1)
        plt.subplots_adjust(top=0.8)
        plt.errorbar(self.decyr[idx], self.N[idx], yerr=self.SN[idx], ecolor='black',
                elinewidth=0.2, capsize=1, capthick=0.5, fmt='o', ms=3, mfc='r', mec='black', mew=0)
        plt.ylabel('North (mm)')

        # East component
        plt.subplot(3,1,2)
        plt.errorbar(self.decyr[idx], self.E[idx], yerr=self.SE[idx], ecolor='black',
                elinewidth=0.2, capsize=1, capthick=0.5, fmt='o', ms=3, mfc='g', mec='black', mew=0)
        plt.ylabel('East (mm)')

        # Vertical component
        plt.subplot(3,1,3)
        plt.errorbar(self.decyr[idx], self.U[idx], yerr=self.SU[idx], ecolor='black',
                elinewidth=0.2, capsize=1, capthick=0.5, fmt='o', ms=3, mfc='b', mec='black', mew=0)
        plt.ylabel('Up (mm)')
        plt.xlabel('Time (year)')
        if len(time_range) == 2:
            plt.xlim(time_range)
        
        plt.suptitle("Time Series of Site Position drawn by zhao at Institute of Seismology\n\nStation: "+self.site+"\n\n %10.3fN %10.3fE %6.2f(m)\n\n %d Daily solution (%7.2f-%7.2f)" %(self.lat, self.lon, self.hei, len(idx), min(self.decyr[idx]), max(self.decyr[idx])), fontsize=15)
        # Adjust the space
        plt.subplots_adjust(hspace=0.4, wspace=0.6)

        if show: plt.show()

        # show the figure
        plt.savefig(self.site+".jpg", format='jpg', dpi=600)



class earthquake(object):
    '''
    Earthquake is a class representing an seismic event.
    '''
    code     = None
    location = []
    epoch    = []
    decyr    = None
    distance = 0.0
    
    def __init__(self, code, location, epoch, distance):
        '''
        Constructor.
        
        Input:
            code      = earthquake ID using 2 char
            location  = [lat, lon, dep] in degree and in km
            epoch     = earthquake time [year, month, day, hour, mininute]
            distance  = distance away from the epicenter in km            
        '''
        self.code     = code
        self.location = location
        self.epoch    = epoch
        self.distance = distance
        jd            = gpstime.ymdhms_to_jd(self.epoch, 0)
        self.decyr    = gpstime.jd_to_decyrs(jd)

        
    def fun_costep(self, t, amp):
        '''
        return offset function
        
        Input:
            t    = a list/array of decimal year
            amp  = amplitude of offset/step
        '''
        return amp * np.heaviside((t-self.decyr), 0)


class eqcatalog(object):
    '''
    eacatalog is a class representing a list of earthquakes.
    '''

    eqlist = []
    def __init__(self, eqfile):
        '''
        Constructor.
        
        Input:
            eqfile  = eq_rename file in GAMIT/GLOBK format
        '''
        
        # check the file exist.
        if os.path.isfile(eqfile) == False:
            logging.info(' The input file %s does not exist!' %(eqfile))
            sys.exit()

        # read eq_rename file
        with open(eqfile) as fid:
            for line in fid:
                if line[0] == ' ':
                    if line.find('eq_def') != -1 or line.find('EQ_DEF') != -1:
                        code = line.split()[1]
                        lat  = float(line.split()[2])
                        lon  = float(line.split()[3])
                        distance = float(line.split()[4])
                        dep  = float(line.split()[5])
                        year = int(line.split()[6])
                        mon  = int(line.split()[7])
                        day  = int(line.split()[8])
                        hour = int(line.split()[9])
                        minu = int(line.split()[10])
                        
                        # eq is an instance of class earthquake
                        eq   = earthquake(code, [lat, lon, dep], 
                                      [year, mon, day, hour, minu], distance )
                        self.eqlist.append(eq)
        return 


    def getEQ(self, code):
        '''
        Return an earthquake instance using an earthquake code

        Input:
            code = 2 upper case char representing an earthquake
        '''
        for i in range(len(self.eqlist)):
            if self.eqlist[i].code == code:
                return self.eqlist[i]
                

class eqPost(object):
    '''
    eqPost is a class representing postseismic term for an earthquake
    '''
    def __init__(self, event, method, mintau, maxtau):
        '''
        Constructor

        Input:
            event  = an instance of class earthquake
            method = 3 char [LOG/EXP]
            mintau = lower bound of tau value should always greater than 0
            maxtau = upper bound of tau value
        '''
        self.eq        = event
        self.method    = method
        self.mintau    = mintau
        self.maxtau    = maxtau



class eqPostList(object):
    '''
    eqPostList is a class representing a list of eqPost 
    '''
    eqpostlist = []
    def __init__(self, eqfile, eqlist):
        '''
        Construvtor

        Input:
            eqfile  = eq_rename.eq
            eqlist  = an instance of class eqcatalog
        '''
        if os.path.isfile(eqfile) == False:
            logging.fatal(' The input file %s does not exist!' %(eqfile))
            sys.exit()
        with open(eqfile) as fid:
            for line in fid:
                if line[0] == " ":
                    if line.find('eq_log') != -1:
                        code   = line.split()[1]
                        mintau = line.split()[2]
                        maxtau = line.split()[3]
                        event  = eqlist.getEQ(code)
                        eqpost = eqPost(event, "LOG", mintau, maxtau)
                        self.eqpostlist.append(eqpost)
                    if line.find('eq_exp') != -1:
                        code   = line.split()[1]
                        mintau = line.split()[2]
                        maxtau = line.split()[3]
                        event  = eqlist.getEQ(code)
                        eqpost = eqPost(event, "LOG", mintau, maxtau)
                        self.eqpostlist.append(eqpost)
        
        
class offset(object):
    '''
    offset is a class representing an break/step in time series.
    '''
    
    def __init__(self, epoch, site):
        '''
        Constructor.
        Mod by Zhao Bin, Jan. 11, 2019. Fix bug of assigning decyr
        
        Input:
            epoch  = break epoch [year, month, day, hour, minute]
            site   = site ID
        '''
        self.epoch = epoch
        self.site  = site
        jd         = gpstime.ymdhms_to_jd(self.epoch, 0)
        self.decyr = gpstime.jd_to_decyrs(jd)

        
    def fun_offset(self, t, amp):
        '''
        return offset function
        
        Input:
            t    = a list/array of decimal year
            amp  = amplitude of offset/step
        '''
        return amp * np.heaviside(t-self.decyr(), 0)


class breakcatalog(object):
    '''
    breakcatalog is a class representing a list instance of class offset
    '''

    breaklist = []
    
    def __init__(self, breakfile):
        '''
        Constructor.
        
        Input:
            breakfile = eq_rename in GAMIT/GLOBK format
        '''
        
        # check the file exist
        if os.path.isfile(breakfile) == False:
            logging.fatal('The file %s does not exist!' %(breakfile))
            sys.exit()
            
        # open and read the file
        with open(breakfile) as fid:
            for line in fid:
                if line[0] != '#':
                    if line.find('break') != -1 or line.find('BREAK') != -1:
                        site = line.split()[1]
                        year = int(line.split()[2])
                        mon  = int(line.split()[3])
                        day  = int(line.split()[4])
                        hour = int(line.split()[5])
                        minu = int(line.split()[6])
                        brk  = offset([year, mon, day, hour, minu], site)
                        self.breaklist.append(brk)


class correction(object):
    '''
    correction is a class representing correction for velocity, breaks, seasonal terms
    '''

    def __init__(self, velfile='', offsetfile='', periodfile=''):
        '''
        Constractor.
        Mod by Zhao Bin, Feb. 17, 2019. Fix bug when reading SITE from correction file.
        Mod by Zhao Bin, Mar.  9, 2020. Decode site name

        Input:
            velfile    = file name of velocity
                         Lon, Lat, E, N, Se, Sn, Cne, Site, U, Su
            offsetfile = file name of offset/break
                         E, N, U, Site, decyr
            periodfile = file name of seasonal terms
                         EAsin, EAcos, ESsin, EScos, NAsin, NAcos, NSsin, NScos, UAsin, UAcos, 
                         USsin, UScos, Site
        '''
        if os.path.isfile(velfile) == False:
            logging.info(" Input velfile %s do not exist!" %(velfile))
            self.correct_velo   = False
        else:
            self.correct_velo   = True
            self.veldata    = np.genfromtxt(velfile, comments='#')
            self.velsite    = np.genfromtxt(velfile, comments='#', usecols=[7], dtype='4S')
            self.velsite    = np.array([i.decode() for i in self.velsite])

        if os.path.isfile(offsetfile) == False:
            logging.info(" Input offsetfile %s do not exist!" %(offsetfile))
            self.correct_offset = False
        else:
            self.correct_offset = True
            self.offsetdata = np.genfromtxt(offsetfile, comments='#', usecols=[0,1,2])
            self.offsetsite = np.genfromtxt(offsetfile, comments='#', usecols=[3], dtype='str')
            self.offsetyear = np.genfromtxt(offsetfile, comments='#', usecols=[4])
#           self.offsetsite = np.array([i.decode() for i in self.offsetsite])

        if os.path.isfile(periodfile) == False:
            logging.info(" Input periodfile %s do not exist!" %(periodfile))
            self.correct_period = False
        else:
            self.correct_period = True
            self.perioddata = np.genfromtxt(periodfile, comments='#')
            self.periodsite = np.genfromtxt(periodfile, comments='#', usecols=[12], dtype='4S')
            self.periodsite = np.array([i.decode() for i in self.periodsite])


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
    
    def __init__(self, site, lon, lat, t, obs, sigma, param_dict, component='N', time_range=[-np.inf, np.inf]):
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
                                len(np.where(abs(self.correct.offsetdata[:,4]-brk.decyr)<0.01)[0]):
                            indx1 = set(np.where(self.correct.offsetsite==self.site)[0])
                            indx2 = set(np.where(abs(self.correct.offsetdata[:,4]-brk.decyr)<0.01)[0])
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
        model = self.ifun(mt, *self.parm)
        plt.plot(self.t, self.obs, 'ro', ms=2)
        plt.plot(self.t, model, color='b')
        
        plt.title(self.site)
        plt.xlabel('Time (year)')
        plt.ylabel('Displacement (mm)')
        plt.show()

    def get_mod(self, time_span=[]):
        '''
        Output modeled time series

        Mod by Zhao Bin, Jan 31, 2019. Adding parameter time_span
        '''
        if len(time_span) == 2:
            mt = np.arange(min(time_span), max(time_span), 1/365.25)
        else:
            mt = np.arange(min(self.t), max(self.t), 1/365.25)
        mfun = self.full_filter(mt)
        m    = mfun(mt, *self.param)
        return mt, m

    def get_correct(self, mod_dict, time_span=[], eqcode=None):
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
        if len(time_span) == 2:
            idx    = np.where(np.logical_and(self.t>time_span[0], self.t<time_span[1]))[0]
            if len(idx) > 0:
                obs_correct = self.ifun(self.t[idx], *param)
        else:
            obs_correct = self.ifun(self.t, *param)
        if len(time_span) == 2:
            mt = np.arange(min(time_span), max(time_span), 1/365.25)
        else:
            mt = np.arange(min(self.t), max(self.t), 1/365.25)
        mfun        = self.full_filter(mt)
        mod_correct = mfun(mt, *param)

        return obs_correct, mod_correct

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
            logging.warning("No correction model for {1:s}".format(nrun.site))
            return
        ndisp = nmod_correct[-1]-nmod_correct[0]

    # east component
    et, em = erun.get_mod()
    idx    = np.where(np.logical_and(et>time_span[0], et<time_span[1]))[0]
    if len(idx) > 0:
        eobs_correct, emod_correct = erun.get_correct(mod_dict, time_span=time_span, eqcode=eqcode)
        if len(emod_correct) == 0: 
            logging.warning("No correction model for {1:s}".format(erun.site))
            return
        edisp = emod_correct[-1]-emod_correct[0]

    # vertical component
    ut, um = urun.get_mod()
    idx    = np.where(np.logical_and(ut>time_span[0], ut<time_span[1]))[0]
    if len(idx) > 0:
        uobs_correct, umod_correct = urun.get_correct(mod_dict, time_span=time_span, eqcode=eqcode)
        if len(emod_correct) == 0: 
            return
            logging.warning("No correction model for {1:s}".format(nrun.site))
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
        outfile = '{}_{}.neu'.format(nrun.site, eqcode)
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
        erun   = an instance of class tsfitting for east component
        urun   = an instance of class tsfitting for up component
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






if __name__ == '__main__':
    import time, getpass
    print("# Created by %s on %s" %(getpass.getuser(), time.asctime()))
    text = "from PyTsfit import *\nimport glob, os\n\neqfile     = ''\nvelfile    = ''\noffsetfile = ''\nperiodfile = ''\n"
    print("%s" %(text))
    text = "eq         = eqcatalg(eqfile)\nbk         = breakcatalog(eqfile)\neqp        = eqPostList(eqfile, eq)"
    print("%s" %(text))
    text = "cor        = correction(velfile, offsetfile, periodfile\n"
    print("%s" %(text))
    text = "plot_dict  = {'detrend':True, 'debreak':True, 'deeqoffset':False, 'depost':False, 'deseason':False}"
    print("%s" %(text))
    text = "mod_dict   = {'detrend':True, 'debreak':True, 'deeqoffset':False, 'depost':False, 'deseason':False}"
    print("%s" %(text))
    text = "param_dict = {\n             'constant'   : True,\n             'linear'     : True,"
    print("%s" %(text))
    text = "             'eqlist'     : eq.eqlist,\n             'eqpostlist' : [],\n             'brklist'    : bk.breaklist"
    print("%s" %(text))
    text = "             'ANN'        : False,\n             'SANN'       : False\n             'correct'    : cor}\n"
    print("%s" %(text))
    print("poslist  = glob.glob('*.pos')")
    print("fid      = open('velo.gmtvec', 'a')")
    print("for posfile in poslist:")
    print("    data = posData(posfile)")
    print("    param_dict['eqpostlist'] = eqp.eqpostlist")
    print("    nrun = tsfitting(data.site, data.lon, data.lat, data.decyr, data.N, data.SN, param_dict, 'N', [2000,2020])")
    print("    nparam = nrun.doFitting()\n")
    print("    param_dict['eqpostlist'] = eqp.eqpostlist")
    print("    erun = tsfitting(data.site, data.lon, data.lat, data.decyr, data.E, data.SE, param_dict, 'E', [2000,2020])")
    print("    eparam = erun.doFitting()\n")
    print("    param_dict['eqpostlist'] = []")
    print("    urun = tsfitting(data.site, data.lon, data.lat, data.decyr, data.U, data.SU, param_dict, 'U', [2000,2020])")
    print("    uparam = erun.doFitting()\n")
    print("    plot_obs_mod(nrun, erun, urun, nparam, eparam, uparam, plot_dict)")
    print("    output_postseismic_disp(nrun, erun, urun, [2015.315, 2017.315], fid=fid)")
    print("    output_postseismic_velo(nrun, erun, urun, [2015.315, 2017.315], fid=fid)")
    print("    output_velo(nrun, erun, urun, fid=fid)")
    print("    output_eqoffset(nrun, erun, urun)")
    print("    output_break(nrun, erun, urun)")
