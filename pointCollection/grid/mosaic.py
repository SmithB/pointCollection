#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mosaic.py
Routines for creating a weighted mosaic from a series of tiles
Written by Tyler Sutterley (03/2020)

UPDATE HISTORY:
    Updated 03/2020: check number of dimensions of z if only a single band
    Written 03/2020
"""

import os
import h5py
import numpy as np
import scipy.ndimage
from .data import data

class mosaic(data):
    def __init__(self, **kwargs):
        #self.x=None
        #self.y=None
        #self.t=None
        #self.z=None
        super().__init__(**kwargs)
        self.mask=None
        self.weight=None
        self.extent=[np.inf,-np.inf,np.inf,-np.inf]
        self.dimensions=[None,None,None]
        self.spacing=[None,None]
        self.fill_value=np.nan

    def update_spacing(self, temp):
        """
        update the step size of mosaic
        """
        self.spacing = (temp.x[1] - temp.x[0], temp.y[1] - temp.y[0])
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
        #if (np.ndim(temp.z) == 3):
        if hasattr(temp,'t') and hasattr(temp.t, 'size') and temp.t.size > 0:
            self.dimensions[2]=temp.t.size
            #ny,nx,self.dimensions[2] = np.shape(temp.z)
            self.t=temp.t.copy()
        else:
            self.dimensions[2] = 1
        # calculate y dimensions with new extents
        self.dimensions[0] = np.int((self.extent[3] - self.extent[2])/self.spacing[1]) + 1
        # calculate x dimensions with new extents
        self.dimensions[1] = np.int((self.extent[1] - self.extent[0])/self.spacing[0]) + 1
        # calculate x and y arrays
        self.x = np.linspace(self.extent[0],self.extent[1],self.dimensions[1])
        self.y = np.linspace(self.extent[2],self.extent[3],self.dimensions[0])
        return self

    def image_coordinates(self, temp):
        """
        get the image coordinates
        """
        iy = np.array((temp.y[:,None]-self.extent[2])/self.spacing[1],dtype=np.int)
        ix = np.array((temp.x[None,:]-self.extent[0])/self.spacing[0],dtype=np.int)
        return (iy,ix)

    def weights(self, pad=0, feather=0, apply=False):
        """
        Create a weight matrix for a given grid
        Apply the weights if specified
        """
        # find dimensions of matrix
        sh = getattr(self, self.fields[0]).shape
        if len(sh)==3:
            nband=sh[2]
        else:
            nband=1
        ny=sh[0]
        nx=sh[1]

        # allocate for weights matrix
        self.weight = np.ones((ny,nx), dtype=np.float)
        gridx,gridy = np.meshgrid(self.x,self.y)
        # pad the weight matrix
        if pad:
            indy,indx = np.nonzero((gridx < (self.x[0] + pad)) |
                (gridx > (self.x[-1] - pad)) |
                (gridy < (self.y[0] + pad)) |
                (gridy > (self.y[-1] - pad)))
            self.weight[indy,indx] = 0.0
        # feathering the weight matrix
        if feather:
            # use a gaussian filter to create smoothed weighting function
            temp = np.ones((ny,nx), dtype=np.float)
            indy,indx = np.nonzero((gridx < (self.x[0] + pad + feather/2)) |
                (gridx > (self.x[-1] - pad - feather/2)) |
                (gridy < (self.y[0] + pad + feather/2)) |
                (gridy > (self.y[-1] - pad - feather/2)))
            temp[indy,indx] = 0.0
            sigma = 0.25*feather/(self.x[1]-self.x[0])
            gauss = scipy.ndimage.gaussian_filter(temp, sigma,
                mode='constant', cval=0)
            # only use points within feather
            indy,indx = np.nonzero((gridx >= (self.x[0] + pad)) &
                (gridx <= (self.x[-1] - pad)) &
                (gridy >= (self.y[0] + pad)) &
                (gridy <= (self.y[-1] - pad)) &
                ((gridx < (self.x[0] + pad + feather)) |
                (gridx > (self.x[-1] - pad - feather)) |
                (gridy < (self.y[0] + pad + feather)) |
                (gridy > (self.y[-1] - pad - feather))))
            self.weight[indy,indx] = gauss[indy,indx]
        # if applying the weights to the original z data
        if apply:
            for field in self.fields:
                if nband > 1:
                    for band in range(nband):
                        getattr(self, field)[:,:,band] *= self.weight
                else:
                    getattr(self, field)[:,:] *= self.weight
        return self

    # def to_h5(self, fileOut, fields, replace=True):
    #     """
    #     write a mosaic object to an hdf5 file
    #     """
    #     # if overwriting the HDF5 file or presently non-existent
    #     if replace or not os.path.isfile(fileOut):
    #         if os.path.isfile(fileOut):
    #             os.remove(fileOut)
    #         fileID=h5py.File(fileOut,'w')
    #     else:
    #         fileID=h5py.File(fileOut,'r+')
    #     # write dimensions to HDF5
    #     h5 = {}
    #     dims = [field for field in fields if field in ('x','y','t')]
    #     for field in dims:
    #         h5[field] = fileID.create_dataset(field, data=getattr(self,field),
    #             compression="gzip")
    #     # write variables to HDF5
    #     for field in sorted(set(fields) - set(dims)):
    #         data = getattr(self,field)
    #         h5[field] = fileID.create_dataset(field, data=data,
    #             fillvalue=self.fill_value, compression="gzip")
    #         # attach dimensions
    #         h5[field].dims[0].label='y'
    #         h5[field].dims[0].attach_scale(h5['y'])
    #         h5[field].dims[1].label='x'
    #         h5[field].dims[1].attach_scale(h5['x'])
    #         if (np.ndim(data) == 3) and ('t' in dims):
    #             h5[field].dims[2].label='t'
    #             h5[field].dims[2].attach_scale(h5['t'])
    #     # close the HDF5 file
    #     fileID.close()
