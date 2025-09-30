#from PyTsfit.tsfitting import tsfitting
from PyTsfit import *
from sklearn import linear_model


import glob
import sys
import numpy as np

sitefile   = 'cmonoc.camp'
sitelist   = np.genfromtxt(sitefile, dtype=str)

for i in range(len(sitelist)):
    posfiles = glob.glob('../pos/'+sitelist[i]+'*.pos')

    for j in range(len(posfiles)):
        posfile = posfiles[j]
        data = posData(posfile)
        print(posfile)
        if len(data.decyr) == 0:
            continue

        ransac = linear_model.RANSACRegressor()


        # North

        year   = data.decyr
        year   = year.reshape((len(year), 1))
        ransac.fit(year, data.N)
        vn     = ransac.estimator_.coef_[0]
        inlier_mask  = ransac.inlier_mask_
        outlier_mask = np.logical_not(inlier_mask)
        if outlier_mask.any() == True:
            print(posfile)
            print(outlier_mask)
