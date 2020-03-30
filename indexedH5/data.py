#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 17:26:25 2019

@author: ben
"""

import numpy as np
import pointCollection as pc
import os
import h5py

class data(pc.data):
    def __init__(self, bin_W=[1.e4, 1.e4], **kwargs):
        self.bin_W=bin_W
        super().__init__(**kwargs)

    def to_file(self, D, out_file, time_field='time', append=True, ind_fields=['x','y','time']):
        y_bin_function=np.round(D.y/self.bin_W[0])
        x_bin_function=np.round(D.x/self.bin_W[1])
        x_scale=np.nanmax(x_bin_function)-np.nanmin(x_bin_function)
        t=getattr(D, time_field)
        t_scale=np.nanmax(t)-np.nanmin(t)

        xy_bin_function=(y_bin_function-np.nanmin(y_bin_function))*x_scale+(x_bin_function-np.nanmin(x_bin_function))
        xyt_bin_function= xy_bin_function + (t-np.nanmin(t))/t_scale
        ind=np.argsort(xyt_bin_function)

        bin_dict={}
        xy_bin_fn_sort=xy_bin_function[ind]
        fn_delta=np.concatenate([[-1], np.flatnonzero(np.diff(xy_bin_fn_sort)), [xy_bin_fn_sort.size]])
        for ii in range(len(fn_delta)-1):
            this_ind=ind[(fn_delta[ii]+1):(fn_delta[ii+1]+1)]
            bin_dict[(x_bin_function[this_ind[0]], y_bin_function[this_ind[0]])]=this_ind
        key_arr=np.array([key for key in bin_dict.keys()])
        key_order=np.argsort(key_arr[:,1]-np.min(key_arr[:,1])*x_scale+(key_arr[:,0]-np.min(key_arr[:,0])))
        key_arr=key_arr[key_order,:]

        for key in key_arr:
            this_group='%dE_%dN' % tuple(key*self.bin_W)
            if append and os.path.isfile(out_file):
                D.copy_subset(bin_dict[tuple(key)]).append_to_h5(out_file, \
                             ind_fields=ind_fields, group=this_group)
            else:
                D.copy_subset(bin_dict[tuple(key)]).to_h5(out_file, \
                             replace=False, group=this_group)

    def read(self, xy_bin, fields=['x','y','time'], index_range=[[-1],[-1]]):
        if isinstance(fields, dict):
            field_list=[]
            for key in fields:
                field_list += fields[key]
        else:
            field_list=fields.copy()
        out_data={field:list() for field in field_list}
        with h5py.File(self.filename,'r') as h5f:
            blank_fields=list()
            if isinstance(xy_bin, np.ndarray):
                xy_bin=[xy_bin[:,0], xy_bin[:,1]]

            if index_range[0][0]>=0:
                # All the geo bins are together.  Use the index_range variable to read
                for field in field_list:
                    if field not in h5f:
                        blank_fields.append(field)
                        continue
                # make sure the index range gets iterated over properly.
                if len(index_range[0])<2:
                    if len(h5f[field].shape)>1:
                        out_data[field].append(np.array(h5f[field][0,int(index_range[0]):int(index_range[1])]))
                    else:
                        out_data[field].append(np.array(h5f[field][int(index_range[0]):int(index_range[1])]))
                else:
                    for i0_i1 in zip(index_range[0], index_range[1]):
                        if len(h5f[field].shape)>1:
                            out_data[field].append(np.array(h5f[field][0,i0_i1[0]:i0_i1[1]]))
                        else:
                            out_data[field].append(np.array(h5f[field][i0_i1[0]:i0_i1[1]]))
            else:
                # this is a file with distinct bins, each with its own set of datasets
                for xy in zip(xy_bin[0], xy_bin[1]):
                    bin_name='%dE_%dN' % xy
                    for field in field_list:
                        if field in h5f:
                            if bin_name in h5f[field]:
                                out_data[field].append(np.array(h5f[field][bin_name]).squeeze())
                        elif bin_name in h5f:
                            if field in h5f[bin_name]:
                                out_data[field].append(np.array(h5f[bin_name][field]).squeeze())
                        else:
                            blank_fields.append(field)
            for field in field_list:
                if isinstance(out_data[field], list):
                    if len(out_data[field])>1 or isinstance(out_data[field], list):
                        try:
                            temp=list()
                            for item in out_data[field]:
                                if item.size > 0 and item.ndim > 0:
                                    temp.append(item)
                                elif item.size==1 and item.ndim ==0:
                                    temp.append(np.array([item]))
                            if len(temp)>1:
                                out_data[field]=np.concatenate(temp)
                            elif len(temp)==0:
                                out_data[field]=np.zeros(0)
                            elif len(temp)==1:
                                out_data[field]=temp[0]
                        except ValueError as e:
                            print("ValueError in read_indexed_h5_file, continuing")
                            print(e)
                else:
                    out_data[field]=np.array(out_data[field])
        return pc.data(fields=field_list).from_dict(out_data)

