#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 14:02:01 2020

@author: ben
"""

import numpy as np
import pointCollection as pc
from .bin_rows import bin_rows

def apply_bin_fn(D_pt, res, fn=None, field='z'):
    
    pt_dict=bin_rows(np.c_[np.round(D_pt.x/res)*res, np.round(D_pt.y/res)*res])
    keys=list(pt_dict.keys())
    
    z=np.array([fn(D_pt, pt_dict[key]) for key in keys])
    N=np.array([len(pt_dict[key]) for key in keys])
    keys=np.array(keys)
    return pc.data(filename=D_pt.filename)\
        .from_dict({'x':keys[:,0], 'y':keys[:,1], field:z,'count':N})
    

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