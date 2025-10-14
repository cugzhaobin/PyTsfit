#!/usr/bin/env python
# Written by Zhao Bin, Institute of Seismology, CEA. Aug 2, 2017

import numpy as np
import matplotlib.pyplot as plt
import os, sys

#timespan = [2018,2022]
nedfiles = [sys.argv[i] for i in range(1, len(sys.argv))]
for nedf in nedfiles:
    dat  = np.genfromtxt(nedf)
    plt.subplot(3,1,1)
    plt.scatter(dat[:,0], dat[:,1], s=2)
#   plt.plot(dat[:,0], dat[:,1], linewidth=2)
    plt.ylabel('N (mm)')
#   plt.xlim(timespan)

    plt.subplot(3,1,2)
    plt.scatter(dat[:,0], dat[:,2], s=2)
#   plt.plot(dat[:,0], dat[:,2], linewidth=2)
    plt.ylabel('E (mm)')
#   plt.xlim(timespan)

    plt.subplot(3,1,3)
    plt.scatter(dat[:,0], dat[:,3], s=2)
#   plt.plot(dat[:,0], dat[:,3], linewidth=2)
    plt.ylabel('U (mm)')
    plt.xlabel('Time (year)')
#   plt.xlim(timespan)

#plt.savefig(nedfiles[0][0:4]+".png", format='png')
#outf = "{}.png".format(nedfiles[0][0:4]) 
#os.popen("eog "+outf)
plt.show()
