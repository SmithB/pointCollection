#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  6 21:11:44 2025

@author: ben
"""

import numpy as np
import pointCollection as pc
import scipy.ndimage as snd
def blend(inputs, dst=None, group=None, field=None, erode=[], feather=None):
    
    temp=dst.copy_meta()
    temp.assign(zw=np.zeros(temp.shape), w=np.zeros(temp.shape))
    dx=dst.x[1]-dst.x[0]

    if len(erode) > 0:
        N_erode = np.ceil(np.array(erode)/dx)
    else:
        N_erode=0

    if np.isscalar(N_erode) or len(N_erode)==1:
        N_erode=np.zeros(len(inputs))+N_erode

    if feather is not None:
        N_feather = np.ceil(feather/dx)

    #loop over inputs
    for count, jj in enumerate(inputs):
        if isinstance(jj, str):
            this = pc.grid.data().from_file(jj, group=group, field_dict={'z':field})
        else:
            this = jj
        this_z = this.interp(dst.x, dst.y, gridded=True)
        mask=np.isfinite(this_z)
        if N_erode[count] > 0:
            mask=snd.binary_erosion(mask, N_erode[count])
            
        if feather is not None:
            wt = snd.distance_transform_edt(mask)
            # check this!
            wt=0.5+0.5*np.cos(np.pi*(1-np.minimum(wt, N_feather)/N_feather))
        else:
            wt=mask.astype(float)
        temp.zw += np.nan_to_num(this_z)*wt
        temp.w += wt
    dst.assign(z=np.zeros(dst.shape)+np.nan)
    nz=temp.w>0
    dst.z[nz] = temp.zw[nz] / temp.w[nz]
    return dst
