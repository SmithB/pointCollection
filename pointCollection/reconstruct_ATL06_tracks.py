#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 11 10:31:24 2020

@author: ben
"""
import pointCollection as pc
import numpy as np


def reconstruct_ATL06_tracks(D, x0=None, y0=None, W=None):
    '''
    sort a collection of data (D) into a list of data organized by track
    '''
    if ('cycle' in D.fields):
        _, bin_dict=pc.unique_by_rows(np.c_[D.cycle, D.rgt, D.BP, D.LR], return_dict=True)
    else:
        _, bin_dict=pc.unique_by_rows(np.c_[D.cycle_number, D.rgt, D.BP, D.LR], return_dict=True)
    D0=[]
    for key in bin_dict:
        if "delta_time" in D.fields:
            ind=np.argsort(D.delta_time[bin_dict[key]])
        else:
            ind=np.argsort(D.time[bin_dict[key]])
        this_D=D[bin_dict[key][ind]]
        if W is not None:
                this_D.index((np.abs(this_D.x-x0)<W/2) &(np.abs(this_D.y-y0)<W/2) )
        D0.append(this_D)   
    return D0
