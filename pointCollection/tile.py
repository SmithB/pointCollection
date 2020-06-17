#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 16:46:30 2019

@author: ben
"""
import numpy as np
import pointCollection as pc
import h5py
import os

class tile(object):
    '''
    A tile is a datafile containing a collection of points from a list of other
    files.  It is written as a pointCollection.indexedH5 file, and contains
    a list of the files that contributed to it, as well as a source_file_num
    dataset that points to each source file
    '''

    def __init__(self, D=None, xy0=None, bin_W=[1.e4, 1.e4], tile_W=1.e5, SRS_proj4=None, time_field=None, z_field=None):
        self.bin_W=bin_W
        self.tile_W=tile_W
        self.SRS_proj4=SRS_proj4
        self.D=D
        self.xy0=xy0
        if time_field is None:
            self.time_field = self.__time_field__()
        if z_field is None:
            self.z_field = self.__z_field__()

    def __default_field_dict__(self):
        return {'None':['x','y','z','time']}

    def __time_field__(self):
        return 'time'

    def __z_field__(self):
        return 'z'

    def from_geoIndex(self, GI_file=None, field_dict=None):

        if field_dict is None:
            field_dict=self.__default_field_dict__()
        dxb, dyb = np.meshgrid(np.arange(-self.tile_W/2, self.tile_W/2+self.bin_W[0], self.bin_W[0]),
                           np.arange(-self.tile_W/2, self.tile_W/2+self.bin_W[1], self.bin_W[1]))
        dxb=dxb.ravel()
        dyb=dyb.ravel()

        fields=[]
        for group in field_dict:
            for ds in field_dict[group]:
                fields.append(ds)

        gI=pc.geoIndex().from_file(GI_file, read_file=False)
        self.D=gI.query_xy((self.xy0[0]+dxb, self.xy0[1]+dyb), fields=field_dict)
        if self.D is not None:
            for Di in self.D:
                if not hasattr(Di,'x'):
                    Di.get_xy(self.SRS_proj4)
                Di.index((np.abs(Di.x-self.xy0[0])<self.tile_W/2) & \
                            (np.abs(Di.y-self.xy0[1])<self.tile_W/2))
        return self

    def write(self, out_dir, fields=None, append=True, ind_fields=['x','y','time']):
        if self.D is None:
            return
        out_file=out_dir+'/E%d_N%d.h5' %(self.xy0[0], self.xy0[1])
        n_source_files=0
        if fields is None:
            field_dict=self.__default_field_dict__()
            fields=[]
            for key in field_dict.keys():
                fields += field_dict[key]
        if os.path.isfile(out_file):
            with h5py.File(out_file,'r') as h5f:
                if 'source_files' in h5f:
                    grp=h5f['source_files']
                    n_source_files=len(h5f['source_files'].attrs.keys())
        file_dict={}

        for file_num, Di in enumerate(self.D):
            if not hasattr(Di,'x'):
                Di.get_xy(self.SRS_proj4)
            Di.assign({'source_file_num':np.zeros_like(Di.x, dtype=int)+file_num+n_source_files})

            Di.ravel_fields()
            Di.index(np.isfinite(getattr(Di, self.z_field)))
            file_dict[file_num]=Di.filename
        out_fields=list(set(fields+['x','y','source_file_num']))
        D_all=pc.data(fields=out_fields).from_list(self.D)
        if D_all.size==0:
            return

        pc.indexedH5.data(filename=out_file, bin_W=self.bin_W).to_file(D_all, \
                         out_file, time_field=self.time_field, append=append,\
                         ind_fields=ind_fields)
        with h5py.File(out_file,'r+') as h5f:
            if 'source_files' in h5f:
                grp=h5f['source_files']
            else:
                grp=h5f.create_group("source_files")
            for key in file_dict:
                grp.attrs['file_%d' % key] = file_dict[key]

