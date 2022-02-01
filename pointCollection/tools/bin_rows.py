#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 13:47:30 2020

@author: ben
"""
import numpy as np

def bin_rows(x):
    """
    determine the unique rows in an array
    
    inputs:
        x:  array of values.  Some function should be applied so that x has
            a limited number of discrete values
    outputs:
        bin_dict: dictionary whose keys are the unique values of x (unsorted) 
        and whose values are vectors of row indices containing those values
    """    
    ind=np.zeros(x.shape[0])
    # bin the data by the values in each column.  From left to right, the importance
    # of the value to the sorting decreases by a factor the number of distinct 
    # values in each column
    scale=1.
    if len(x.shape)==1:
        x.shape=[x.shape[0], 1]
    for col in range(x.shape[1]):        
        z, ii=np.unique(x[:,col].astype(np.float64), return_inverse=True)
        scale /= (np.max(ii).astype(float)+1.)
        ind += ii * scale     
    u_ii, index, inverse=np.unique(ind, return_index=True, return_inverse=True)

    bin_dict={}
    inv_arg_ind=np.argsort(inverse)
    inv_sort=inverse[inv_arg_ind]
    ind_delta=np.concatenate([[-1], np.where(np.diff(inv_sort))[0], [inv_sort.size]])
    for ii in range(len(ind_delta)-1):
        this_ind=inv_arg_ind[(ind_delta[ii]+1):(ind_delta[ii+1]+1)]
        bin_dict[tuple(x[this_ind[0], :])]=this_ind
    return bin_dict