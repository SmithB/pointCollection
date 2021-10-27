#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 14:02:01 2020

@author: ben
"""

import numpy as np
import pointCollection as pc
from .bin_rows import bin_rows

def apply_bin_fn(D_pt, res, fn=None, fields=['z'], xy0=[0, 0]):

    pt_dict=bin_rows(np.c_[np.round((D_pt.x-xy0[0])/res)*res+xy0[0], np.round((D_pt.y-xy0[1])/res)*res+xy0[1]])
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
        XR=[np.min(x.ravel()), np.max(x.ravel())]
        YR=[np.min(y.ravel()), np.max(y.ravel())]
        print(XR)
        xg=np.arange(XR[0], XR[1], res)
        yg=np.arange(YR[0], YR[1], res)
        zg=np.zeros((len(yg), len(xg)))+background
        grid=pc.grid.data().from_dict({'x':xg, \
                                     'y':yg, 
                                     'z':zg})
    else:
        XR=grid.extent[0:2]
        YR=grid.extent[2:4]
        zg=grid.z.copy()
    c=((x-XR[0])/res).astype(int)
    r=((y-YR[0])/res).astype(int)
    ii = (c>0) & (c<zg.shape[1])
    ii &=  (r>0) & (r<zg.shape[0])
    zg[tuple(np.c_[r[ii].ravel(), c[ii].ravel()].T)]=\
        getattr(D_pt, field)[ii].ravel()
    return grid