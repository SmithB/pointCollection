#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 12:07:54 2020

@author: ben
"""

import os
import re
import contextlib
import numpy as np
import pointCollection as pc

# ATL11_ttttrr_c0c1_rrr_vv.h5 -- tttt is the reference ground track (RGT),
# equal to /ancillary_data/start_rgt, so rgt can be read from the filename
# with no file I/O at all.
_ATL11_FILENAME_RE = re.compile(r'^ATL11_(\d{4})\d{2}_\d{4}_\d{3}_\d{2}\.h5$')


class data(pc.data):
    def __init__(self, pair=2, **kwargs):
        self.pair=pair
        self.pair_name = f'pt{int(self.pair)}'
        self.cycle_number=np.array([])
        super().__init__(**kwargs)


    def __default_field_dict__(self, field_weight='light'):
        if field_weight=='light':
            field_dict={None:['latitude','longitude','h_corr',\
                                   'h_corr_sigma', 'h_corr_sigma_systematic',\
                                   'delta_time','quality_summary', 'ref_pt'], \
                    'ref_surf': ['dem_h', 'x_atc','fit_quality']}
        return self.__convert_field_dict__(field_dict)

    def __convert_field_dict__(self, field_dict):

        temp={}
        for key in field_dict:
            if key is None:
                temp[self.pair_name]=field_dict[key]
            elif key in ['__calc_internal__']:
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
            # The only field that should be tiled by row is 'cycle_number'.
            # If others come up they can be added to this list
            if field in ['cycle_number']:
                if temp.ndim==1:
                        temp.shape=(1,cols)
                setattr(self, field, np.tile(temp, [rows,1]))
            elif temp.size==rows:
                if temp.ndim==1:
                        temp.shape=(rows,1)
                # BUGFIX: don't tile empty fields that already have the right number of columns
                if temp.size==0 and np.max(temp.shape)==cols:
                    continue
                setattr(self, field, np.tile(temp, [1, cols]))
        self.__update_size_and_shape__()
        
        return self

    def __internal_field_calc__(self, field_dict, fs=None, h5_f=None):
        if 'rgt' in field_dict['__calc_internal__']:
            self.__update_size_and_shape__()
            m = _ATL11_FILENAME_RE.match(os.path.basename(self.filename)) if self.filename else None
            if m is not None:
                # no file access needed at all -- the RGT is encoded in the filename
                self.rgt=int(m.group(1))+np.zeros(self.shape)
            elif h5_f is not None:
                self.rgt=h5_f['/ancillary_data/start_rgt'][0]+np.zeros(self.shape)
            else:
                with pc.io_utils.open_h5(self.filename, fs=fs) as h5f:
                    self.rgt=h5f['/ancillary_data/start_rgt'][0]+np.zeros(self.shape)

    def from_h5(self, filename, pair=None, field_weight='light', tile_fields=True, fs=None, h5_f=None, **kwargs):
        if pair is not None:
            self.pair=pair
            self.pair_name=f'pt{int(pair)}'
            self.field_dict=self.__default_field_dict__(field_weight=field_weight)
        _ctx = contextlib.nullcontext(h5_f) if h5_f is not None else pc.io_utils.open_h5(filename, fs=fs)
        with _ctx as h5f:
            cycle_number = np.array(h5f[self.pair_name]['cycle_number'])

            if 'field_dict' in kwargs and kwargs['field_dict'] is not None:
                kwargs['field_dict']=self.__convert_field_dict__(kwargs['field_dict'].copy())

            super().from_h5(filename, h5_f=h5f, fs=fs, **kwargs)
        self.cycle_number=cycle_number
        self.columns=len(self.cycle_number)
        self.shape=(self.latitude.size, self.columns)
        if tile_fields:
            self.fields += ['cycle_number']
            self.__tile_fields__()
        return self
