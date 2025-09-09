#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 08:36:04 2024

@author: ben
"""
import scipy
import numpy as np


def smooth_corrected(z, mask, w_smooth, **filter_args):

    filter_defaults= {"mode":"constant", "cval":0}
    if filter_args is None:
        filter_args={}
    else:
        for arg, val in filter_args.items():
            filter_defaults[arg]=val

    mask1=scipy.ndimage.gaussian_filter(mask.astype(float), w_smooth, **filter_defaults)
    ztemp=np.nan_to_num(z)
    ztemp[mask==0]=0.0

    zs=scipy.ndimage.gaussian_filter(ztemp, w_smooth, **filter_defaults)
    zs[mask1>0]=zs[mask1>0]/mask1[mask1>0]
    return zs, mask1

def fill_edges(z, w_smooth=1, ndv=-9999, dim=None, **filter_args):

    z1=z.copy()
    if dim is None:
        # assume 2-d input
        temp=z1
        mask=(temp != ndv) & (np.isfinite(temp))
        temp[mask==0]=np.nan
        temp_s, mask1=smooth_corrected(temp, mask, w_smooth, **filter_args)
        temp_s[mask1==0]=np.nan
        temp[mask==0]=temp_s[mask==0]
        z1=temp
    elif dim==2:
        for ii in range(z.shape[2]):
            temp=z1[:,:,ii]
            mask=(temp != ndv) & (np.isfinite(temp))
            temp[mask==0]=np.nan
            temp_s, mask1=smooth_corrected(temp, mask, w_smooth, **filter_args)
            temp_s[mask1==0]=np.nan
            temp[mask==0]=temp_s[mask==0]
            z1[:,:,ii]=temp
    elif dim==0:
        for ii in range(z.shape[1]):
            temp=z1[ii,:,:]
            mask=(temp != ndv) & (np.isfinite(temp))
            temp[mask==0]=np.nan
            temp_s, mask1=smooth_corrected(temp, mask, w_smooth, **filter_args)
            temp_s[mask1==0]=np.nan
            temp[mask==0]=temp_s[mask==0]
            z1[ii,:,:]=temp
    elif dim ==1:
        raise IndexError('dimension 1 not implemented')
    return z1
