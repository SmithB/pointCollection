#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  1 08:39:06 2025

@author: ben
"""

import numpy as np
import os
import pointCollection as pc
import json
import re
import glob

class tilingSchema(object):
    def __init__(self, tile_spacing=1.e5, tol=None,
                 mapping_function_name='round',
                 mapping_function=np.round,
                 EPSG=None,
                 coords=['x','y'],
                 scale=1000,
                 format_str='E%d_N%d',
                 extension='.h5',
                 bin_size=1.e4,
                 directory=None):
        self.tile_spacing = tile_spacing
        self.bin_size=bin_size
        if mapping_function is not None:
            self.mapping_function = mapping_function
            mapping_function_name = self.mapping_function.__name__
        self.mapping_function_name=mapping_function_name
        self.extension = extension
        self.EPSG = EPSG
        self.coords = coords
        self.format_str = format_str
        self.format_re  = re.compile(self.format_str.replace(r'%d',r'(.*)')+self.extension)
        self.scale = scale
        self.data_format = 'indexedH5'
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
                      'scale', 'format_str', 'extension','directory','bin_size']:
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
                if 'tiling_schema' in fh:
                    group='tiling_schema'
                else:
                    group='/'
                for key, val in fh[group].items():
                    scheme_dict[key] = val
        for key, val in scheme_dict.items():
            if hasattr(self, key):
                setattr(self, key, val)
        if self.directory is None:
            self.directory = os.path.dirname(scheme_file)
        return self

    # TBD: implement latlon keyword
    def tile_xy(self, all_tiles=True,
                unique=True,
                return_dict=False,
                xy=None,
                data=None, tol=None):

        if tol is None:
            tol=self.bin_size/2

        if self.mapping_function is None:
            self.set_mapping_function()

        if xy is None:
            xy = [data.x.ravel(), data.y.ravel()]

        if np.isscalar(xy[0]):
            xy = [*map(np.atleast_1d, xy)]

        if len(xy[0]) > 0 and not isinstance(xy, np.ndarray):
            xy = [*map(np.array,xy)]

        if return_dict:
            # return one tile xy for each point:
            tile_xys = self.mapping_function( np.c_[xy[0], xy[1]] / self.tile_spacing ) * self.tile_spacing
            _, tile_dict = pc.unique_by_rows(tile_xys, return_dict=True)
            return tile_dict

        # return the unique tile centers that could
        # contribute to the points specified by xy0
        if all_tiles and self.mapping_function_name=='round':
            # need to check for xys that are on boundaries.  For those that are, add
            # another point that is just on the other side of the boundary
            for dim, other_dim in zip([0, 1], [1, 0]):
                for sgn in [-1, 1]:
                    ctrs = np.round(xy[dim]/self.tile_spacing)*self.tile_spacing
                    delta =  xy[dim] - ctrs
                    # check for points at the upper end of this bin
                    bdry_ind = np.flatnonzero(sgn * delta >= 0.5*self.tile_spacing - tol)
                    xy[dim] = np.append(xy[dim], ctrs[bdry_ind] + sgn*(self.tile_spacing/2 + tol), axis=0)
                    xy[other_dim] = np.append(xy[other_dim], xy[other_dim][bdry_ind])
        tile_xys = self.mapping_function( np.c_[xy[0], xy[1]] / self.tile_spacing ) * self.tile_spacing
        if unique:
            return np.unique(tile_xys, axis=0)
        else:
            return tile_xys

    def tile_filename(self, xy_t):
        return os.path.join(self.directory, self.format_str % tuple([xy_t[0]/self.scale, xy_t[1]/self.scale]))+self.extension

    def filenames_for_xy(self, xy0):
        if np.isscalar(xy0[0]):
            xy0=[np.array([xy0[0]]).ravel(), np.array([xy0[1]]).ravel()]
        if not isinstance(xy0[0], np.ndarray):
            xy0=[*map(np.array, xy0)]
        tile_filenames = []
        for xyt in self.tile_xy(xy=xy0, unique=True, all_tiles=True):
            tile_filenames.append(self.tile_filename(xyt))
        return tile_filenames

    def filenames_for_box(self, xyr, resolution=1.e4):
        xg, yg = np.meshgrid( np.arange(xyr[0][0], xyr[0][1] + resolution * 1.01, resolution),
                              np.arange(xyr[1][0], xyr[1][1] + resolution * 1.01, resolution) )
        return self.filenames_for_xy([xg.ravel(), yg.ravel()])

    def tile_bounds(self, xy = [0.,0.]):
        if self.mapping_function==np.round:
            offset = [0,0]
        elif self.mapping_function == np.floor:
            offset = [self.tile_spacing/2, self.tile_spacing/2]
        xyT = self.tile_xy(xy=xy)[0]
        return [xy_i + off_i + np.array([-1, 1])*self.tile_spacing/2 for xy_i, off_i in zip(xyT, offset)]

    def tile_boundary(self, xy = [0., 0.]):
        bds = self.tile_bounds(xy)
        return (bds[0][[0, 0, 1, 1, 0]], bds[1][[0, 1, 1, 0, 0]])

    def write_tiles(self, D, bin_size=None, replace=True):
        tile_dict = self.tile_xy(data=D, return_dict=True)
        for xy0, ii in tile_dict.items():
            out_file = self.tile_filename(xy0)
            if self.data_format == 'h5':
                D[ii].to_h5(out_file, replace=True)
            elif self.data_format == 'indexed_h5':
                pc.indexedH5.data( bin_W = (bin_size, bin_size) ).to_file(D[ii], out_file, replace=replace)

    def file_xy(self, filenames=None):
        if filenames is None:
            filenames=glob.glob(os.path.join(self.directory,'*.h5'))
        xy=[]
        for filename in filenames:
            m = self.format_re.search(os.path.basename(filename))
            if m is not None:
                xy.append(np.array([*map(float, m.groups())])*self.scale)
        return xy
