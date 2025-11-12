#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  6 21:11:44 2025

@author: ben
"""

import numpy as np
import pointCollection as pc
import scipy
def blend(inputs, dst=None, inplace=False, group=None, field='z', erode=[], feather=None):

    if dst is None:
        dst=inputs[0].copy_meta()

    temp=dst.copy_meta()
    temp.assign(zw=np.zeros(temp.shape), w=np.zeros(temp.shape))
    dx=dst.x[1]-dst.x[0]

    if erode is not None:
        N_erode = np.ceil(np.array(erode)/dx)
    else:
        N_erode=0

    if np.isscalar(N_erode) or len(N_erode)==1:
        N_erode=np.zeros(len(inputs))+N_erode

    print(N_erode)

    if feather is not None:
        N_feather = np.ceil(feather/dx)

    #loop over inputs
    for count, jj in enumerate(inputs):
        if isinstance(jj, str):
            this = pc.grid.data().from_file(jj, group=group, field_dict={'z':field})
        else:
            this = jj
        this_z = this.interp(dst.x, dst.y, field=field, gridded=True)
        mask=np.isfinite(this_z)

        if N_erode[count] > 0:
            kernel=np.ones(np.array([N_erode[count], N_erode[count]]).astype(int))
            mask=scipy.ndimage.binary_erosion(mask, kernel, border_value=1)

        if feather is not None:
            wt = scipy.ndimage.distance_transform_edt(np.pad(mask,1))[1:-1,:][:, 1:-1]
            # check this!
            wt=0.5+0.5*np.cos(np.pi*(1-np.minimum(wt, N_feather)/N_feather))
        else:
            wt=mask.astype(float)
        temp.zw += np.nan_to_num(this_z)*wt
        temp.w += wt
    temp.assign(z=np.zeros(dst.shape)+np.nan)
    nz=temp.w>0
    temp.z[nz] = temp.zw[nz] / temp.w[nz]
    dst.assign({field:temp.z})
    if inplace:
        return
    else:
        return dst
