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
        self.spacing=spacing
        self.fill_value=np.nan

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
        if hasattr(temp,'t') and hasattr(temp.t, 'size') and temp.t.size > 0:
            self.dimensions[2]=temp.t.size
            self.t=temp.t.copy()
        else:
            self.dimensions[2] = 1
        # calculate y dimensions with new extents
        self.dimensions[0] = np.int64((self.extent[3] - self.extent[2])/self.spacing[1]) + 1
        # calculate x dimensions with new extents
        self.dimensions[1] = np.int64((self.extent[1] - self.extent[0])/self.spacing[0]) + 1
        # calculate x and y arrays
        self.x = self.extent[0] + self.spacing[0]*np.arange(self.dimensions[1])
        self.y = self.extent[2] + self.spacing[1]*np.arange(self.dimensions[0])
        return self

    def image_coordinates(self, temp):
        """
        get the image coordinates
        """
        iy = np.array((temp.y[:,None]-self.extent[2])/self.spacing[1],dtype=np.int64)
        ix = np.array((temp.x[None,:]-self.extent[0])/self.spacing[0],dtype=np.int64)
        return (iy,ix)

    def raised_cosine_weights(self, pad, feather):
        """
        Create smoothed weighting function using a raised cosine function
        """
        weights=[]
        for dim, xy in zip([0, 1], [self.x, self.y]):
            xy0 = np.mean(xy)
            W = xy[-1]-xy[0]
            dist = np.abs(xy-xy0)
            wt = np.zeros_like(dist)
            i_feather = (dist >= W/2 - pad - feather) & ( dist <= W/2 -pad )
            wt_feather = 0.5 + 0.5 * np.cos( -np.pi*(dist[i_feather] - (W/2 - pad - feather))/feather)
            wt[ i_feather ] = wt_feather
            wt[ dist <= W/2 - pad - feather ] = 1
            wt[ dist >= W/2 - pad] = 0
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

    def weights(self, pad=0, feather=0, apply=False, mode='raised cosine'):
        """
        Create a weight matrix for a given grid
        Apply the weights if specified
        """
        # find dimensions of matrix
        sh = getattr(self, self.fields[0]).shape
        if len(sh)==3:
            ny, nx, nband = sh
        else:
            nband =1
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
        # if applying the weights to the original z data
        if apply:
            for field in self.fields:
                if nband > 1:
                    for band in range(nband):
                        getattr(self, field)[:,:,band] *= self.weight
                else:
                    getattr(self, field)[:,:] *= self.weight
        return self
