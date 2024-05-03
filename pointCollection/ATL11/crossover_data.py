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
                 'latitude', 'longitude', 'ref_pt', 'rgt','tide_ocean','dac']}
        return {self.pair_name+'/'+key:field_dict[key] for key in field_dict}


    def from_h5(self, filename, pair=None, D_at=None, crossover_fields=None):
        '''
        Read crossover data from a file.  

        inputs:
            filename: the file to read
            pair: the pair to read from the file
            D_at (optional): along-track data already read from the file.  Only
                the reference points within D_at will be read.
        outputs:
            self: a data structure with fields of shape n_points x n_cycles x 2
              The first dimension gives different crosover points
              The second dimension has one entry for each cycle in the dataset
              The third dimension separates reference-point data (index 0)
                 and crossing-track data (index 1)
        '''
        if D_at is not None:
            # make the pair match the pair for D_at
            pair=D_at.pair
            
        if pair is None:
            pair=self.pair
        else:
            self.pair=pair
 
        self.pair_name = f'pt{int(self.pair)}'

        if D_at is None:
            D_at=pc.ATL11.data().from_h5(filename, pair=pair)
            if len(D_at.ref_pt)==0:
                return None
            index_range=None
        else:
            # copy D_at to avoid side effects
            D_at=D_at.copy()
            if len(D_at.ref_pt)==0:
                return None
            ref_pt_range = [np.nanmin(D_at.ref_pt), np.nanmax(D_at.ref_pt)]
            temp=pc.data().from_h5(filename, field_dict={self.pair_name+'/crossing_track_data':['ref_pt']})
            ind=np.flatnonzero((temp.ref_pt >= ref_pt_range[0]) & (temp.ref_pt <=ref_pt_range[1]))
            if len(ind)==0:
                return None
                #index_range=[None, None]
            else:
                index_range=[np.min(ind), np.max(ind)]
        if crossover_fields is None:
            field_dict = self.__default_XO_field_dict__()
        else:
            field_dict = {self.pair_name+'/crossing_track_data':crossover_fields}
        D_xo=pc.data().from_h5(filename, field_dict=field_dict, index_range=index_range)
        D_xo.index(np.isfinite(D_xo.ref_pt) & np.isfinite(D_xo.cycle_number) )
        with h5py.File(filename,'r') as h5f:
            rgt=int(h5f[self.pair_name].attrs['ReferenceGroundTrack'])

        D_at.assign({'rgt':np.zeros_like(D_at.delta_time)+rgt})
        n_cycles=D_at.cycle_number[-1, -1]

        u_pt_xo=np.unique(D_xo.ref_pt)
        D_at.index(np.in1d(D_at.ref_pt[:,0], u_pt_xo))
        theshape=(u_pt_xo.size, n_cycles,2)
        all_fields=set(D_xo.fields+D_at.fields)
        self.assign({field:np.zeros(theshape)+np.NaN for field in all_fields})

        row=np.searchsorted( u_pt_xo, D_at.ref_pt.ravel())
        col=D_at.cycle_number.ravel().astype(int)-1    
        self_ind0=np.ravel_multi_index((row, col, np.zeros_like(row, dtype=int)), theshape)
               
        for field in D_at.fields:
            getattr(self, field).flat[self_ind0]=getattr(D_at, field).ravel()

        row=np.searchsorted(u_pt_xo, D_xo.ref_pt)
        col=D_xo.cycle_number.astype(int)-1
        self_ind1=np.ravel_multi_index((row, col, 1+np.zeros_like(row, dtype=int)), theshape)
        for field in self.fields:
            if field in D_xo.fields:
                getattr(self, field).flat[self_ind1]=getattr(D_xo, field)

        self.__update_size_and_shape__()
        return self

