#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mosaic.py
Routines for creating a weighted mosaic from a series of tiles

UPDATE HISTORY:
    Updated 06/2023: calculate x and y arrays using np.arange and spacing
    updated 03/2021: change scheme for calculating weights, raised cosine as default
    Updated 03/2020: check number of dimensions of z if only a single band
    Written 03/2020
"""

import numpy as np
from .data import data
import pointCollection as pc

class mosaic(data):
    def __init__(self, spacing=[None,None], **kwargs):
        #self.x=None
        #self.y=None
        #self.t=None
        #self.z=None
        super().__init__(**kwargs)
        self.invalid=None
        self.weight=None
        self.extent=[np.inf,-np.inf,np.inf,-np.inf]
        self.dimensions=[None,None,None]
        self.field_dims={}
        self.spacing=spacing
        self.tile_weight=None
        self.fill_value=np.nan
        self.normalized=True
        self.use_time=False
        self.fields=[]

    def from_grid(self, source, copy=False):
        """make a mosaic object from a pc.data.grid object"""

        for field in ['x','y','projection','filename','extent','time', 't', 't_axis']:
            if hasattr(source, field):
                setattr(self, field, getattr(source, field))
        for field in source.fields:
            if copy:
                self.assign({field:getattr(source, field).copy()})
            else:
                self.assign({field:getattr(source, field)})
        self.__update_size_and_shape__()
        self.__update_extent__()
        return self

    def update_spacing(self, temp):
        """
        update the step size of mosaic
        """
        # try automatically getting spacing of tile
        try:
            self.spacing = (temp.x[1] - temp.x[0], temp.y[1] - temp.y[0])
        except:
            pass
        return self

    def update_bounds(self, temp):
        """
        update the bounds of mosaic
        """
        if (temp.extent[0] < self.extent[0]):
            self.extent[0] = np.copy(temp.extent[0])
        if (temp.extent[1] > self.extent[1]):
            self.extent[1] = np.copy(temp.extent[1])
        if (temp.extent[2] < self.extent[2]):
            self.extent[2] = np.copy(temp.extent[2])
        if (temp.extent[3] > self.extent[3]):
            self.extent[3] = np.copy(temp.extent[3])
        return self

    def update_dimensions(self, temp):
        """
        update the dimensions of the mosaic with new extents
        """
        # get number of bands
        t_name=None
        for this_t_name in['t','time']:
            if hasattr(temp, this_t_name):
                t_attr=getattr(temp, this_t_name)
                if hasattr(t_attr, '__len__') and len(t_attr) > 0:
                    self.dimensions[2]=len(t_attr)
                    setattr(self, this_t_name, t_attr.copy())
                    t_name=this_t_name
        if t_name is None:
            self.dimensions[2] = 1
        # calculate y dimensions with new extents
        self.dimensions[0] = np.int64((self.extent[3] - self.extent[2])/self.spacing[1]) + 1
        # calculate x dimensions with new extents
        self.dimensions[1] = np.int64((self.extent[1] - self.extent[0])/self.spacing[0]) + 1
        # calculate x and y arrays
        self.x = self.extent[0] + self.spacing[0]*np.arange(self.dimensions[1])
        self.y = self.extent[2] + self.spacing[1]*np.arange(self.dimensions[0])
        return self

    def image_coordinates(self, temp, validate=False):
        """
        get the image coordinates
        """
        iy = np.rint((temp.y[:,None]-self.extent[2])/self.spacing[1]).astype(np.int64)
        ix = np.rint((temp.x[None,:]-self.extent[0])/self.spacing[0]).astype(np.int64)

        if validate:
            iy1 = np.flatnonzero((iy >= 0) & (iy < self.shape[0]))[:, None]
            iy0 = iy[iy1.ravel(),:]
            ix1 = np.flatnonzero((ix >= 0) & (ix < self.shape[1]))[None,:]
            ix0 = ix[:,ix1.ravel()]
            return(iy0, ix0, iy1, ix1)
        else:
            return (iy,ix)

    def setup_bounds_from_list(self, in_list,
                group=None,
                fields=None,
                bounds=None,
                spacing=None):

        for item in in_list.copy():
            if isinstance(item, str):
                # read tile grid from file
                try:
                    temp=pc.grid.data().from_file(item, group=group, meta_only=True)
                    if bounds is not None:
                        temp=temp.cropped(*bounds)
                    if temp is not None and (len(temp.x)>0) and (len(temp.y) > 0):
                        #update the mosaic bounds to include this tile
                        self.update_bounds(temp)
                        if self.spacing[0] is None:
                            self.update_spacing(temp)
                    else:
                        in_list.remove(item)
                except Exception:
                    print(f"failed to read group {group} "+ str(item))
                    in_list.remove(item)
            else:
                if bounds is not None:
                    item=item.cropped(*bounds)
                if item is not None and (len(item.x)>0) and (len(item.y) > 0):
                    self.update_bounds(item)
                    if self.spacing[0] is None:
                        self.update_spacing(item)
                    temp=item
                else:
                    in_list.remove(item)

        self.update_dimensions(temp)
        self.__update_extent__()
        self.__update_size_and_shape__()


    def setup_fields(self, item, group=None, fields=None):
        '''
        Set up fields based on an input data item
        '''

        # create output mosaic, weights, and mask
        # read data grid from the first tile HDF5, use it to set the field dimensions

        if isinstance(item, str):
            prototype=pc.grid.mosaic().from_file(item, group=group, fields=fields)
        else:
            prototype=item
        if len(prototype.fields) == 0:
            message = f"pointCollection.grid.mosaic.py: did not find fields {fields} in file {prototype.filename}"
            return message
        # if time is not specified, squeeze extra dimenstions out of inputs
        self.use_time=False
        for field in ['time','t']:
            if hasattr(prototype, field) and getattr(prototype, field) is not None:
                self.use_time=True
        if not self.use_time:
            for field in prototype.fields:
                setattr(prototype, field, np.squeeze(getattr(prototype, field)))
        if fields is None:
            fields=prototype.fields
        these_fields=[field for field in fields if field in prototype.fields]
        self.field_dims={field:getattr(prototype, field).ndim for field in these_fields}
        for field in these_fields:
            self.assign({field:np.zeros(self.dimensions[0:self.field_dims[field]])})
        self.invalid = np.ones(self.dimensions,dtype=bool)
        self.__update_size_and_shape__()


    def replace(self, item, group=None, fields=None):
        """
        Overwrite a section of the mosaic with an input file or mosaic
        """

        if fields is None:
            fields=self.fields.copy()
        # read data grid from HDF5
        if isinstance(item, str):
            temp=pc.grid.mosaic().from_file(item, group=group, fields=fields)
        else:
            temp=item
        if not temp.overlaps(self):
            return

        these_fields=[field for field in fields if field in temp.fields]
        # get the image coordinates of the input file
        iy0, ix0, iy1, ix1 = self.image_coordinates(temp, validate=True)
        for field in these_fields:
            try:
                if self.field_dims[field]==3:
                    field_data=np.atleast_3d(getattr(temp, field))
                    for band in range(self.dimensions[2]):
                        getattr(self, field)[iy0,ix0,band] = field_data[iy1,ix1,band]
                    self.invalid[iy0,ix0] = False
                else:
                    field_data=getattr(temp, field)
                    getattr(self, field)[iy0,ix0] = field_data[iy1, ix1]
                    self.invalid[iy0, ix0] = False
            except IndexError as e:
                thestr = f"problem with field {field}"
                if isinstance(item, str):
                    thestr += "in group {group} in file {item}"
                else:
                    thestr += "in item {in_list.index(item))}"
                print(thestr)
                raise(e)


    def add(self, item, fields, group=None, use_time=False,
            pad=0, feather=0):
        """
        Add all bands from an item to a mosaic.

        This method propagates invalid values from the input to the output, and
        adds all bands at once

        Parameters
        ----------
        item : str or mosaic
            Item to be added to the mosaic.
        fields : iterable
            Fields to be added.
        group : str optional
            Group from which to read fields, if item is a file  The default is None.
        use_time : bool, optional
            If true, preserve singleton dimensions in inputs. The default is False.
        pad : float, optional
            pad the weights by this distance. The default is 0.
        feather : float optional
            feather weights by this distance. The default is 0.

        Returns
        -------
        None.

        """

        if fields is None:
            fields=self.fields.copy()

        self.normalized=False
        # read data grid from file
        if isinstance(item, str):
            temp=pc.grid.mosaic().from_file(item, group=group, fields=fields)
        else:
            if isinstance(item, self.__class__):
                temp=item
            else:
                temp=pc.grid.mosaic().from_grid(item)
        if not temp.overlaps(self):
            return
        temp.update_spacing(temp)
        if not use_time:
            for field in temp.fields:
                setattr(temp, field, np.squeeze(getattr(temp, field)))
        these_fields=[field for field in fields if field in temp.fields]
        # copy weights for tile
        if self.tile_weight is not None:
            temp.weight=self.tile_weight.copy()
        else:
            temp.weights(pad=pad, feather=feather)
            self.tile_weights=temp.weight.copy()
        # get the image coordinates of the input file
        iy0, ix0, iy1, ix1 = self.image_coordinates(temp, validate=True)
        for field in these_fields:
            try:
                if self.field_dims[field]==3:
                    field_data=np.atleast_3d(getattr(temp, field))
                    bands=range(self.dimensions[2])
                    for band in bands:
                        getattr(self, field)[iy0,ix0,band] += field_data[iy1,ix1,band]*temp.weight[iy1, ix1]
                    self.invalid[iy0,ix0] = False
                else:
                    field_data=getattr(temp, field).copy()
                    getattr(self, field)[iy0, ix0] += field_data*temp.weight[iy1,ix1]
                    self.invalid[iy0,ix0] = False
            except (IndexError, ValueError) as e:
                thestr = f"problem with field {field}"
                if isinstance(item, str):
                    thestr += "in group {group} in file {item}"
                else:
                    thestr += "in item {in_list.index(item))}"
                print(thestr)
                raise(e)
        # add weights to total weight matrix
        self.weight[iy0,ix0] += temp.weight[iy1,ix1]


    def add_to_band(self, item, fields, group=None,
                        pad=0, feather=0, band=None, in_band=None, out_band=None):
        """
        Add an item to a band of a mosaic.

        This method treats invalids in inputs as providing no data, rather
        than propagating the invalid to the output

        Parameters
        ----------
        item : str or mosaic
            Item to be added to the mosaic.
        fields : iterable
            Fields to be added.
        group : str optional
            Group from which to read fields, if item is a file  The default is None.
        use_time : bool, optional
            If true, preserve singleton dimensions in inputs. The default is False.
        pad : float, optional
            pad the weights by this distance. The default is 0.
        feather : float optional
            feather weights by this distance. The default is 0.
        in_band: band to read from input object
        out_band: band in self to write

        Returns
        -------
        None.

        """
        if in_band is None:
            in_band=band

        if out_band is None and in_band is not None:
            out_band = in_band

        if fields is None:
            fields=self.fields.copy()

        self.normalized=False
        # read data grid from file
        if isinstance(item, str):
            if in_band is None:
                temp=pc.grid.mosaic().from_file(item, group=group, fields=fields)
            else:
                temp=pc.grid.mosaic().from_file(item, group=group, fields=fields, bands=[in_band])
        else:
            if isinstance(item, self.__class__):
                if in_band is None:
                    temp=item
                else:
                    temp=item[:,:,in_band]
            else:
                if in_band is None:
                    temp=pc.grid.mosaic().from_grid(item)
                else:
                    temp=pc.grid.mosaic().from_grid(item[:,:,in_band])
        temp.update_spacing(temp)
        if not temp.overlaps(self):
            return

        for field in temp.fields:
            setattr(temp, field, np.squeeze(getattr(temp, field)))
        these_fields=[field for field in fields if field in temp.fields]
        # copy weights for tile
        if self.tile_weight is not None:
            temp.weight=self.tile_weight.copy()
        else:
            temp.weights(pad=pad, feather=feather)
            self.tile_weight=temp.weight.copy()
        # get the image coordinates of the input file
        iy0, ix0, iy1, ix1 = self.image_coordinates(temp, validate=True)
        for field in these_fields:
            try:
                field_data=getattr(temp, field).copy()
                valid_mask = np.isfinite(field_data)
                if not np.all(valid_mask):
                    temp.weight[valid_mask==0]=0
                    field_data[valid_mask==0]=0
                getattr(self, field)[iy0, ix0, out_band] += field_data[iy1, ix1] * temp.weight[iy1,ix1]
                self.invalid[iy0, ix0] = False
            except (IndexError, ValueError) as e:
                thestr = f"problem with field {field}"
                if isinstance(item, str):
                    thestr += "in group {group} in file {item}"
                else:
                    thestr += "in item {in_list.index(item))}"
                print(thestr)
                raise(e)
        # add weights to total weight matrix

        self.weight[iy0,ix0] += temp.weight[iy1,ix1]

    def normalize(self, fields=None, band=None, by_weight=True):
        """
        Normalize the mosaic fields by the sum of the weights

        Returns
        -------
        None.

        """

        # find invalid points:
        if by_weight:
            i_zero= np.flatnonzero( (self.weight == 0) | self.invalid)
            i_nonzero = np.flatnonzero(self.weight)
        else:
            i_zero= np.flatnonzero( self.invalid )

        # find valid points
        for field in self.fields:
            if self.field_dims[field]==3 and band is None:
                band_list = range(self.dimensions[2])
            elif band is not None:
                band_list=[band]
            else:
                band_list=[None]
            for this_band in band_list:
                if this_band is None:
                    if by_weight:
                        getattr(self, field)[:, :, this_band].ravel()[i_nonzero] /= self.weight.ravel()[i_nonzero]
                    getattr(self, field)[:, :, this_band].ravel()[i_zero] = self.fill_value
                else:
                    if by_weight:
                        iy_nz, ix_nz = np.unravel_index(i_nonzero, self.shape[0:2])
                        i_out = np.ravel_multi_index((iy_nz, ix_nz, np.zeros_like(ix_nz)+this_band), self.shape)
                        #temp_out=getattr(self, field)[:, :, this_band]
                        #temp_out.ravel()[i_nonzero] /= self.weight.ravel()[i_nonzero]
                        #getattr(self, field)[:,:,this_band]=temp_out
                        getattr(self, field).ravel()[i_out] /= self.weight.ravel()[i_nonzero]
                    if len(i_zero) > 0:
                        iy_z, ix_z = np.unravel_index(i_zero, self.shape[0:2])
                        i_out = np.ravel_multi_index((iy_z, ix_z, np.zeros_like(ix_z)+this_band), self.shape)
                        getattr(self, field).ravel()[i_out] = self.fill_value

        if self.weight is not None:
            self.weight[:]=0
        self.normalized=True

    def raised_cosine_weights(self, pad, feather):
        """
        Create smoothed weighting function using a raised cosine function
        """
        weights=[]
        for dim, xy, delta in zip([0, 1], [self.y, self.x], self.spacing):
            xy0 = np.mean(xy)
            eps_grid=0.01*delta
            W = xy[-1]-xy[0]
            dist = np.abs(xy-xy0)
            wt = np.zeros_like(dist)
            i_feather = (dist >= W/2 - pad - feather - eps_grid) & ( dist <= W/2 -pad )
            wt_feather = 0.5 + 0.5 * np.sin( -np.pi*(dist[i_feather] - (W/2 - pad - feather/2)) / (feather+2*delta))
            wt[ i_feather ] = wt_feather
            wt[ dist < W/2 - pad - feather - eps_grid] = 1
            wt[ dist > W/2 - pad + eps_grid] = 0
            weights += [wt]
        self.weight *= weights[0][:,None].dot(weights[1][None,:])

    def gaussian_weights(self, pad, feather):
        """
        Create smoothed weighting function using a Gaussian function
        """
        weights=[]
        for dim, xy in zip([0, 1], [self.x, self.y]):
            xy0 = np.mean(xy)
            W = xy[-1]-xy[0]
            dist = np.abs(xy-xy0)
            wt = np.zeros_like(dist)
            i_feather = (dist >= W/2 - pad - feather) & ( dist <= W/2 -pad )
            wt_feather = np.exp(-((xy[i_feather]-xy0)/(feather/2.))**2)
            wt[ i_feather ] = wt_feather
            wt[ dist <= W/2 - pad - feather ] = 1
            wt[ dist >= W/2 - pad] = 0
            weights += [wt]
        self.weight *= weights[0][:,None].dot(weights[1][None,:])

    def pad_edges(self, pad):
        """
        Pad the edges of the weights with zeros
        """
        weights=[]
        for dim, xy in zip([0, 1], [self.x, self.y]):
            xy0 = np.mean(xy)
            W = xy[-1]-xy[0]
            dist = np.abs(xy-xy0)
            wt=np.ones_like(dist)
            wt[ dist >= W/2 - pad] = 0
            weights += [wt]
        self.weight *= weights[0][:,None].dot(weights[1][None,:])

    def weights(self, pad=0, feather=0, mode='raised cosine'):
        """
        Create a weight matrix for a given grid
        """
        # find dimensions of matrix
        sh = getattr(self, self.fields[0]).shape
        if len(sh)==3:
            ny, nx, nband = sh
        else:
            ny, nx = sh
        # allocate for weights matrix
        self.weight = np.ones((ny,nx), dtype=float)
        # feathering the weight matrix
        if feather:
            if mode == 'raised cosine':
                self.raised_cosine_weights(pad, feather)
            elif mode == 'gaussian':
                self.gaussian_weights(pad, feather)
        if pad:
            self.pad_edges(pad)

        return self

    def from_list(self, in_list,
                  bounds=None,
                  fields=None,
                  group='/',
                  pad=0,
                  feather=0,
                  by_band=True,
                  ):
        """
        Generate a mosaic from a list of inputs.

        Inputs can be strings (indicating files) or pointCollection.grid or
        pointCollection.mosaic objects.

        Parameters
        ----------
        in_list : iterable
            Objects to mosaic.
        bounds : iterable of iterables, optional
            specifies bounds of output mosaic, [[xmin, xmax], [ymin, ymax]].
            if None, the bounds will be determined from the inputs
            The default is None.
        fields : iterable of str, optional
            Specifies fields from inputs to be mosaicked.  If not specified, all
            will be mosaicked.  The default is None.
        pad : float, optional
            elements within pad will be removed from inputs before mosacking.
            The default is 0.
        feather : float, optional
            Smooth blending length for inputs. The default is 0.
        group : str, optional
            group in hdf5 or netcdf4 files to read. The default is '/'.

        Returns
        -------
        self
            pointCollection.grid.mosaic object containing the mosaicked data.

        """
        weight = (pad is not None and pad > 0) or (feather is not None and feather>0)

        self.setup_bounds_from_list(in_list, group=group, fields=fields, bounds=bounds)
        message = self.setup_fields(in_list[0], group=group, fields=fields)
        if message is not None:
            return message
        # check if using a weighted summation scheme for calculating mosaic
        if weight:
            if by_band:
                if len(self.shape)>2:
                    band_list=range(self.shape[2])
                else:
                    band_list=[None]
                for band in band_list:
                    self.invalid = np.ones(self.dimensions[0:2],dtype=bool)
                    self.weight = np.zeros((self.dimensions[0],self.dimensions[1]))
                    for item in in_list:
                        self.add_to_band(item, group=group, fields=fields, pad=pad, feather=feather, band=band)
                    self.normalize(band=band)
            else:
                self.invalid = np.ones(self.dimensions[0:2],dtype=bool)
                self.weight = np.zeros((self.dimensions[0],self.dimensions[1]))
                # for each file in the list
                for item in in_list:
                   self.add(item, group=group, fields=fields, pad=pad, feather=feather)
                self.normalize()
        else:
            # overwrite the mosaic with each subsequent tile
            # for each file in the list
            self.invalid = np.ones(self.dimensions[0:2],dtype=bool)
            for item in in_list:
                self.replace(item, group=group, fields=fields)
            self.normalize(by_weight=False)

        return self
