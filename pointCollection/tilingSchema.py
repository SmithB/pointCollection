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
                 mapping_function=None,
                 EPSG=None,
                 coords=['x','y'],
                 scale=1000,
                 format_str='E%d_N%d',
                 format_variables=['x','y'],
                 extension='.h5',
                 bin_size=1.e4,
                 tile_offset = [0,0],
                 directory=None,
                 source=None):
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
        self.format_variables = format_variables
        self.scale = scale
        self.data_format = 'indexedH5'
        self.directory = directory
        self.tile_offset = tile_offset
        self.mapping_function = mapping_function
        # optional dict describing a remote source for the tiles, e.g.
        # {'type': 'EarthAccess', 'short_name': 'ATL11XO', 'daac': 'NSIDC'}.
        # When set, tile_filename() returns bare granule-name strings
        # instead of joining them to a local `directory`, and
        # resolve_files_for_box() resolves those names to real URLs via a
        # single batched earthaccess.search_data(granule_name=[...]) call.
        self.source = source

    def set_mapping_function(self, mapping_function_name=None):

        if mapping_function_name is None:
            mapping_function_name = self.mapping_function_name
        if mapping_function_name == 'round':
            self.mapping_function = np.round
        elif mapping_function_name == 'floor':
            self.mapping_function = np.floor
        else:
            raise NotImplementedError(f'mapping function {mapping_function_name} not understood')
        self.mapping_function_name = self.mapping_function.__name__

    def _scheme_dict(self):

        scheme_dict={}
        for field in ['tile_spacing','mapping_function_name', 'EPSG', 'coords',
                      'tile_offset', 'format_str','format_variables',
                      'scale', 'extension','directory','bin_size','source']:
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
        if self.directory is None and self.source is None:
            self.directory = os.path.dirname(scheme_file)
        return self

    # TBD: implement latlon keyword
    def tile_xy(self, all_tiles=True,
                unique=True,
                return_dict=False,
                xy=None,
                data=None, tol=None):

        # break out the offset attribute to a numpy array
        xy0 = np.array(self.tile_offset).ravel()

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
            tile_xys = xy0 + \
                self.mapping_function(
                    np.c_[xy[0]-self.tile_offset[0], xy[1] - self.tile_offset[1]]
                        / self.tile_spacing ) * self.tile_spacing
            _, tile_dict = pc.unique_by_rows(tile_xys, return_dict=True)
            return tile_dict

        # return the unique tile centers that could
        # contribute to the points specified by xy0
        if all_tiles and self.mapping_function_name=='round':
            # need to check for xys that are on boundaries.  For those that are, add
            # another point that is just on the other side of the boundary
            for dim, other_dim in zip([0, 1], [1, 0]):
                for sgn in [-1, 1]:
                    ctrs = np.round((xy[dim]-xy0[dim])/self.tile_spacing)*self.tile_spacing + xy0[dim]
                    delta =  xy[dim] - ctrs
                    # check for points at the upper end of this bin
                    bdry_ind = np.flatnonzero(sgn * delta >= 0.5*self.tile_spacing - tol)
                    xy[dim] = np.append(xy[dim], ctrs[bdry_ind] + sgn*(self.tile_spacing/2 + tol), axis=0)
                    xy[other_dim] = np.append(xy[other_dim], xy[other_dim][bdry_ind])
        tile_xys = self.mapping_function(
                        (np.c_[xy[0], xy[1]]-xy0) / self.tile_spacing ) * self.tile_spacing + xy0
        if unique:
            return np.unique(tile_xys, axis=0)
        else:
            return tile_xys

    def tile_filename(self, xy_t):
        if set(['xmin','xmax','ymin','ymax']) == set(self.format_variables):
            var_val = {var:val for var, val in zip(['xmin','xmax','ymin','ymax'],
                                                np.concatenate(self.tile_bounds(xy_t)))}
            vals_sorted  = [(var_val[var]/self.scale)
                                for var in self.format_variables]
            name = self.format_str % tuple(vals_sorted) + self.extension
        elif  set(['x','y']) == set(self.format_variables):
            xy0 = [xy_t[0]/self.scale, xy_t[1]/self.scale]
            name = self.format_str % tuple(xy0) + self.extension
        else:
            return None
        if self.source is not None:
            # remote source: return the bare granule/file name to search
            # for, not a local path (there is no `directory` to join)
            return name
        return os.path.join(self.directory, name)

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

    def resolve_files_for_box(self, xyr, fs=None, resolution=1.e4, verbose=False):
        """
        Find the tiles overlapping a box, and resolve each to a location
        that can actually be opened -- a local path, an S3 URI already
        confirmed to exist, or (if self.source specifies a remote source)
        a URL resolved via a search against that source. Tiles that can't
        be found are silently dropped (reported if verbose=True).

        Parameters
        ----------
        xyr : 2-element iterable of 2-element iterables
            [[xmin, xmax], [ymin, ymax]] box bounds.
        fs : s3fs.S3FileSystem, optional
            filesystem to use / reuse for remote existence checks and reads.
            If None and self.source specifies EarthAccess, one is obtained
            via pc.io_utils.get_s3fs() and returned for the caller to reuse.
        resolution : float, optional
            grid spacing used to enumerate candidate tile centers within
            the box (passed to filenames_for_box()).
        verbose : bool, optional
            print a message for each candidate tile that isn't found.

        Returns
        -------
        resolved : dict
            {tile_basename: resolved_location}, for tiles that were found.
        fs : s3fs.S3FileSystem or None
            the filesystem used (for the caller to reuse on subsequent calls).
        """
        candidates = self.filenames_for_box(xyr, resolution=resolution)
        if self.source is not None and self.source.get('type') == 'EarthAccess':
            import earthaccess
            search_kwargs = {k: v for k, v in self.source.items()
                              if k not in ('type', 'daac')}
            earthaccess.login(strategy='netrc')
            granules = earthaccess.search_data(granule_name=candidates, **search_kwargs)
            found = {}
            for g in granules:
                url = g.data_links(access='direct')[0]
                found[os.path.basename(url)] = url
            resolved = {name: found[name] for name in candidates if name in found}
            if verbose:
                for name in candidates:
                    if name not in found:
                        print(f'tilingSchema: {name} not found via earthaccess')
            if fs is None:
                fs = pc.io_utils.get_s3fs(daac=self.source.get('daac', 'NSIDC'))
        else:
            resolved = {}
            for name in candidates:
                if pc.io_utils.path_exists(name, fs=fs):
                    resolved[os.path.basename(name)] = name
                elif verbose:
                    print(f'tilingSchema: {name} not found')
        return resolved, fs

    def tile_bounds(self, xy = [0.,0.]):
        if self.mapping_function is None:
            self.set_mapping_function()
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
                if set(['x','y']) == set(self.format_variables):
                    xy.append(np.array([*map(float, m.groups())])*self.scale)
                elif ['xmin','xmax','ymin','ymax'] == self.format_variables:
                    this_xy =np.array([*map(float, m.groups())])*self.scale
                    xy.append(np.array(np.mean(this_xy[0:2]),
                                       np.mean(this_xy[2:4])))
        return xy

# example:
# ts = tilingSchema(format_str='z0%d_%d_%d_%d', format_variables=['xmin','xmax','ymin','ymax'], tile_spacing=2.e5, EPSG=3031, tile_offset=[1.e5, 1.e5])
