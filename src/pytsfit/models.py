#!/usr/bin/env python
# ----------------------------------------------------------
# Event and constraint models.
#
# earthquake / eqcatalog / eqPost / eqPostList : seismic events
# offset   / breakcatalog                      : non-earthquake steps
# correction                                   : velocity/offset/seasonal priors
#
# Extracted from PyTsfit.py during the module split. No behaviour change.
# ----------------------------------------------------------
import os, sys, logging
import numpy as np
from . import GPSTime as gpstime

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
