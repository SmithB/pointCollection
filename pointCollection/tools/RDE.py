# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 16:31:30 2017

@author: ben
"""

import numpy as np

def RDE(x):
    xs=x.copy()
    xs=np.isfinite(xs)   # this changes xs from values to a boolean
    if np.sum(xs)<2 :
        return np.nan
    ind=np.arange(0.5, np.sum(xs))
    LH=np.interp(np.array([0.16, 0.84])*np.sum(xs), ind, np.sort(x[xs]))
    #print('LH =',LH)
    return (LH[1]-LH[0])/2.  # trying to get some kind of a width of the data ~variance


#import scipy.stats as stats
#def RDE(x):
#    return (stats.scoreatpercentile(x, 84 )-stats.scoreatpercentile(x, 16))/2.