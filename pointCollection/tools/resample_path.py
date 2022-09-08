# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 09:40:49 2019

@author: ben
"""
import numpy as np
def resample_path(x, y, spacing):
    """
    interpolate a path to a different resolution
    
    inputs:
        x, y : path coordinates
        spacing: intended resolution of the result
        
    Calculates the along-track distance for each point in x and y, then interpolates
    the x and y coordinates to an evenly spaced vector of distance values, separated by
    'spacing'
    """     
    
    s0=np.concatenate([[0], np.cumsum(np.diff(np.abs(x+1j*y)))])
    s=np.arange(s0[0], s0[-1]+spacing, spacing)
    if s[-1] < s0[-1]:
        s=np.concatenate([s, [s0[-1]]])
    return np.interp(s, s0, x), np.interp(s, s0, y)