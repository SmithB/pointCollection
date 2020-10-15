#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 14:02:01 2020

@author: ben
"""

import numpy as np
import pointCollection as pc
from .bin_rows import bin_rows

def apply_bin_fn(D_pt, res, fn=None, fields=['z']):

    pt_dict=bin_rows(np.c_[np.round(D_pt.x/res)*res, np.round(D_pt.y/res)*res])
    keys=list(pt_dict.keys())
    xy=np.array(keys)
    result=pc.data(filename=D_pt.filename)\
        .from_dict({'x':xy[:,0].ravel(), 'y':xy[:,1].ravel(), 'count':np.array([len(pt_dict[key]) for key in keys])})
    for field in fields:
        result.assign({field:np.zeros(len(keys), dtype=float)})
    z=np.zeros((len(keys), len(fields)))
    for ii, key in enumerate(keys):
        temp=fn(D_pt, pt_dict[key])
        if len(fields)==1:
            z[ii]=temp
        else:
            z[ii,:]=temp
    for col, field in enumerate(fields):
        result.assign({field:z[:, col].ravel()})
    return result
    

def points_to_grid(D_pt, res, grid=None, field='z', background=np.NaN):
    x=np.round(D_pt.x/res)*res
    y=np.round(D_pt.y/res)*res
    
    if grid is None:
        XR=[np.min(x), np.max(x)]
        YR=[np.min(y), np.max(y)]
        zg=np.zeros((np.diff(YR)/res, np.diff(XR)/res))+background
        grid=pc.grid.data().from_dict({'x':np.arange(XR[0], XR[1], res), \
                                     'y':np.arange(YR[0], YR[1], res), 
                                     'z':zg}, projection=D_pt.projection)
    else:
        XR=grid.extent[0:2]
        YR=grid.extent[2:4]
        zg=grid.z
    c=np.int(x-XR[0])/res
    r=np.int(y-YR[0])/res
    ii = (c>0) & (c<zg.shape[0])
    ii &=  (r>0) & (r<zg.shape[1])
    zg[tuple(np.c_[c[ii].ravel(), r[ii].ravel()].T)]=\
        getattr(D_pt, field)[ii].ravel()
    return grid