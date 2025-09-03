#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  1 08:39:06 2025

@author: ben
"""

import numpy as np
import pointCollection as pc
import os
import json

class tilingSchema(object):
    def __init__(self, tile_spacing=1.e5, mapping_function_name='round',
                 mapping_function=None, EPSG=None, coords=['x','y'],
                 scale=1000, format_str='E%d_N%d', extension='.h5',
                 directory=None):
        self.tile_spacing = tile_spacing
        if mapping_function is not None:
            self.mapping_function = mapping_function
            mapping_function_name = self.mapping_function.__name__
        self.mapping_function_name=mapping_function_name
        self.EPSG = EPSG
        self.coords = coords
        self.format_str = format_str
        self.extension = extension
        self.scale = scale
        self.directory = directory
        self.mapping_function = mapping_function

    def set_mapping_function(self, mapping_function_name=None):

        if mapping_function_name is None:
            mapping_function_name = self.mapping_function_name
        if mapping_function_name == 'round':
            self.mapping_function = np.round
        elif mapping_function_name == 'floor':
            self.mapping_function = np.floor
        self.mapping_function_name = self.mapping_function.__name__

    def _scheme_dict(self):

        scheme_dict={}
        for field in ['tile_spacing','mapping_function_name', 'EPSG', 'coords',
                      'scale', 'format_str', 'extension','directory']:
            try:
                scheme_dict[field] = float(getattr(self, field))
            except (ValueError, TypeError):
                scheme_dict[field] = getattr(self, field)
        return scheme_dict

    def print(self):
        for key, val in self._scheme_dict().items():
            print(f'{key} -> {str(val)}')

    def to_json(self, json_file):

        scheme_dict=self._scheme_dict()
        with open(json_file,'w') as fh:
            json.dump(scheme_dict, fh, indent=2)

    def from_file(self, scheme_file):

        # choose what kind of file this is:
        if scheme_file.endswith('.json'):
            with open(scheme_file,'r') as fh:
                scheme_dict = json.load(fh)
        elif scheme_file.endswith('.h5'):
            import h5py
            scheme_dict={}
            with h5py.File(scheme_file,'r') as fh:
                for key, val in fh.items():
                    scheme_dict[key] = val
        for key, val in scheme_dict.items():
            if hasattr(self, key):
                setattr(self, key, val)
        if self.directory is None:
            self.directory = os.path.dirname(scheme_file)
        return self

    # TBD: implement latlon keyword
    def tile_xy(self, xy=None, data=None, return_dict=False):
        if self.mapping_function is None:
            self.set_mapping_function()
        if xy is None:
            xy = [data.x.ravel(), data.y.ravel()]
        if np.isscalar(xy[0]):
            xy = [xy[0], xy[1]]
        tile_xys = self.mapping_function( np.c_[*xy] / self.tile_spacing ) * self.tile_spacing            
        return pc.unique_by_rows(tile_xys, return_dict=return_dict) 

    def tile_filename(self, xy_t):
        return os.path.join(self.directory, self.format_str % tuple([xy_t[0]/self.scale, xy_t[1]/self.scale]))+self.extension

    def filenames_for_xy(self, xy0):
        if np.isscalar(xy0[0]):
            xy0=[np.array([xy0[0]]).ravel(), np.array([xy0[1]]).ravel()]
        if not isinstance(xy0[0], np.ndarray):
            xy0=[*map(np.array, xy0)]
        return [self.tile_filename(xyt) for xyt in self.tile_xy(xy0)]


    def filenames_for_box(self, xyr, resolution=1.e4):
        xg, yg = np.meshgrid( np.arange(xyr[0][0], xyr[1][1] + resolution * 1.01, resolution),
                              np.arange(xyr[0][0], xyr[1][1] + resolution * 1.01, resolution) )
        return self.filenames_for_xy([xg.ravel(), yg.ravel()])
