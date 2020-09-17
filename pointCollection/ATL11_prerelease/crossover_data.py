#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 13:27:26 2020

@author: ben
"""

import pointCollection as pc
import numpy as np
import h5py

class crossover_data(pc.data):
    def __init__(self, pair=2, **kwargs):
        self.pair=pair
        self.pair_name = f'pt{int(self.pair)}'
        self.columns=2
        super().__init__(**kwargs)


    def __default_XO_field_dict__(self):

        field_dict= {'crossing_track_data':\
                ['along_track_rss', 'atl06_quality_summary', 'cycle_number',
                 'delta_time', 'h_corr', 'h_corr_sigma', 'h_corr_sigma_systematic',
                 'latitude', 'longitude', 'ref_pt', 'rgt']}
        return {self.pair_name+'/'+key:field_dict[key] for key in field_dict}


    def from_h5(self, filename, pair=None):
        if pair is None:
            pair=self.pair
        else:
            self.pair_name = f'pt{int(self.pair)}'
        D_at=pc.ATL11.data().from_h5(filename, pair=pair)
        D_xo=pc.data().from_h5(filename, group=self.pair_name+'/crossing_track_data', field_dict=self.__default_XO_field_dict__())
        D_xo.index(np.isfinite(D_xo.ref_pt) & np.isfinite(D_xo.cycle_number))
        with h5py.File(filename,'r') as h5f:
            rgt=int(h5f[self.pair_name].attrs['ReferenceGroundTrack'])

        D_at.assign({'rgt':np.zeros_like(D_at.delta_time)+rgt})
        n_cycles=D_at.cycle_number[-1, -1]

        u_pt_xo=np.unique(D_xo.ref_pt)
        D_at.index(np.in1d(D_at.ref_pt[:,0], u_pt_xo))

        theshape=(u_pt_xo.size, n_cycles,2)
        self.assign({field:np.zeros(theshape)+np.NaN for field in D_xo.fields})

        row=np.searchsorted( u_pt_xo, D_at.ref_pt.ravel())
        col=D_at.cycle_number.ravel().astype(int)-1
        for field in self.fields:
            if field not in D_at.fields:
                continue
            getattr(self, field).flat[np.ravel_multi_index((row, col, np.zeros_like(row, dtype=int)), theshape)]=getattr(D_at, field).ravel()

        row=np.searchsorted(u_pt_xo, D_xo.ref_pt)
        col=D_xo.cycle_number.astype(int)-1
        for field in self.fields:
            getattr(self, field).flat[np.ravel_multi_index((row, col, 1+np.zeros_like(row, dtype=int)), theshape)]=getattr(D_xo, field)

        self.__update_size_and_shape__()
        return self

