#!/usr/bin/env python
# ----------------------------------------------------------
# Data containers that read a time series file.
#
# neuData : a GAMIT/GLOBK .neu file
# posData : a PBO .pos file
#
# Extracted from PyTsfit.py during the module split. No behaviour change.
# ----------------------------------------------------------
import os, sys, logging
import numpy as np
from . import geotools as gt
from . import GPSTime as gpstime
import matplotlib.pyplot as plt

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
