#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 10:32:16 2022

@author: ben
"""

import numpy as np


def interp_pts_from_grid(D0, Dg, field='z', band=None):
    ''' interpolate points from a pc.grid object

    Parameters
    ----------
    D0 : list pc.data or dict or geodataframe
        object containing x and y coordinates for which the interpolation will be calculated
        Possibilities include:
            -two -field iterable containing x and y numpy arrays
            -pointCollection.data or other container with x and y fields
            -dict or dataframe or other container with x and y entries
    Dg : pointCollection.grid.data
        grid containing the data to be interpolated
    field : str, optional
        field to interpolated. The default is 'z'.
    band : int, optional
        the band in Dg to interpolate. The default is None.

    Raises
    ------
    AttributeError
        Error thrown if attempts to extract x and y from D0 fail

    Returns
    -------
    v : TYPE
        DESCRIPTION.

    '''

    # try to work out what D0 is
    if isinstance(D0, (list, tuple) ):
        x=D0[0].copy()
        y=D0[1].copy()
    else:
        try:
            # maybe D0 is a dict or a dataFrame
            x=D0['x'].copy()
            y=D0['y'].copy()
        except AttributeError:
            # it wasn't a dict or a dataFrame
            try:
                # maybe D0 is a pointCollection object or a namedTuple (or something)
                x=D0.x.copy()
                y=D0.y.copy()
            except AttributeError:
                raise AttributeError('interp_pts_from_grid could not find x and y in the girst argument')
            pass

    # extra-fast interpolator for a regular grid
    ccf=(x-Dg.x[0])/(Dg.x[1]-Dg.x[0])
    rrf=(y-Dg.y[0])/(Dg.x[1]-Dg.x[0])
    good=(ccf >= 0) & (ccf <= Dg.shape[1]-2) &  (rrf >= 0) & (rrf <= Dg.shape[0]-2)
    ccf=ccf[good]
    rrf=rrf[good]

    rr0=np.floor(rrf).astype(int)
    cc0=np.floor(ccf).astype(int)
    rr=np.c_[rr0, rr0+1, rr0, rr0+1]
    cc=np.c_[cc0, cc0, cc0+1, cc0+1]
    v=np.zeros(x.shape)+np.nan
    if band is None:
        z=getattr(Dg, field)

    elif Dg.t_axis==2 or Dg.t_axis is None:
        z=getattr(Dg, field)[:,:,band]
    elif Dg.t_axis==0:
        z=getattr(Dg, field)[band,:,:]

    # sum the the values multiplied by their bilinear weights.  the 'clip' option
    # means that indices past the edge of the grid return the values at the edge
    # of the grid.

    v[good]=np.sum(
            z.ravel()[np.ravel_multi_index([rr.ravel(), cc.ravel()], z.shape)].reshape(rr.shape)*\
                (1-np.abs(np.tile(rrf[:,None], [1, 4])-rr))*(1-np.abs(np.tile(ccf[:, None], [1, 4])-cc)),\
                axis=1, mode='clip')

    return v

