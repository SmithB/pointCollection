#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 12:07:54 2020

@author: ben
"""

import numpy as np
import pointCollection as pc
import os
import h5py


class data(pc.data):
    def __init__(self, pair=2, **kwargs):
        self.pair=pair
        self.pair_name = f'pt{int(self.pair)}'
        self.cycle_number=np.array([])
        super().__init__(**kwargs)


    def __default_field_dict__(self, field_weight='light'):
        if field_weight=='light':
            field_dict={'corrected_h':['latitude','longitude','h_corr',\
                                   'h_corr_sigma', 'h_corr_sigma_systematic',\
                                   'delta_time','quality_summary', 'ref_pt'], \
                    'ref_surf': ['dem_h', 'x_atc']}
        return self.__convert_field_dict__(field_dict)

    def __convert_field_dict__(self, field_dict):

        temp={}
        for key in field_dict:
            if key is None:
                temp[self.pair_name]=field_dict[key]
            if key in ['__calc_internal__']:
                # skip appending the pair name to the __calc_internal__ field
                temp[key]=field_dict[key]
            else:
                temp[self.pair_name+'/'+key]=field_dict[key]
        return temp


    def __tile_fields__(self):
        rows=self.shape[0]
        cols=self.shape[1]
        for field in self.fields:
            temp=getattr(self, field)
            if temp.size==cols:
                if temp.ndim==1:
                        temp.shape=(1,cols)
                setattr(self, field, np.tile(temp, [rows,1]))
            elif temp.size==rows:
                if temp.ndim==1:
                        temp.shape=(rows,1)
                setattr(self, field, np.tile(temp, [1, cols]))
        self.__update_size_and_shape__()
        return self

    def __internal_field_calc__(self, field_dict):
        if 'rgt' in field_dict['__calc_internal__']:
            self.__update_size_and_shape__()
            with h5py.File(self.filename,'r') as h5f:
                self.rgt=h5f['/ancillary_data/start_rgt'][0]+np.zeros(self.shape)

    def from_h5(self, filename, pair=None, field_weight='light', tile_fields=True, **kwargs):
        if pair is not None:
           self.pair_name=f'pt{int(pair)}'
           self.field_dict=self.__default_field_dict__(field_weight=field_weight)
        with h5py.File(filename,'r') as h5f:
            cycle_number = np.array(h5f[self.pair_name]['corrected_h']['cycle_number'])

        if 'field_dict' in kwargs:
            kwargs['field_dict']=self.__convert_field_dict__(kwargs['field_dict'].copy())

        super().from_h5(filename, **kwargs)
        self.cycle_number=cycle_number
        self.columns=len(self.cycle_number)
        self.shape=(self.latitude.size, self.columns)
        if tile_fields:
            self.fields += ['cycle_number']
            self.__tile_fields__()
        return self
