#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 14:13:07 2022

@author: ben
"""
import numpy as np

def seg_dist_cplx(xy, d2_xy, atc_xy, xy0_p, dxy_p, s_p, ds_p, d_max):

    # calculate the potential of each point
    dATC=np.conj(dxy_p)*(xy-xy0_p)/ds_p
    ii = (np.real(dATC) >= -d_max) & (np.real(dATC) <= ds_p + d_max) & \
        (np.abs(np.imag(dATC)) < d_max )
    dATC=dATC[ii]
    d2_temp = np.imag(dATC)**2
    these=np.real(dATC)<0
    d2_temp[these] += np.real(dATC[these])**2
    these=np.real(dATC) > ds_p
    d2_temp[these] += (np.real(dATC[these])-ds_p)**2
    smaller = d2_temp < d2_xy[ii]
    if np.any(smaller):
        ii[ii] = smaller
        #pdb.set_trace()
        d2_xy[ii] = d2_temp[smaller]
        atc_xy[ii] = s_p + dATC[smaller]

def along_track_coords(xy0, xyp, d_max=None, make_grid=False):
    '''
    Calculate the along-track coordinates of a set of points relative to a path

    Parameters
    ----------
    xy0 : iterable
        two coordinate arrays for which the along-track coordinates will be
        calculated.  If 'make_grid' is true, these will be converted to grids
    xyp : iterable
        two coordinate arrays for the path
    d_max : float, optional
        maximum distance from the path for which coordinates will be calculated.
        Points farther than this from the path will be assigned coordinates of
        np.nan
    make_grid : TYPE, optional
        If true, meshgrid will be applied to the coordinates in xy0.
        The default is False.

    Returns
    -------
    at_full : numpy array
        Along-track coordinates for points in xy0.
    xt_full : TYPE
        Across-track coordinates for points in xy0.
    dist : TYPE
        Distance from the path to each point.

    '''
    if make_grid:
        gx, gy=np.meshgrid(xy0[0], xy0[1])
        xy=gx+1j*gy
    else:
        xy=xy0[0]+1j*xy0[1]

    if d_max is not None:
        XYR = [(np.nanmin(xypi)-d_max, np.nanmax(xypi)+d_max) for xypi in xyp]
        ind = np.ones_like(xy, dtype=bool)
        for fn, XYRi in zip((np.real, np.imag), XYR):
            ind &= (fn(xy) > XYRi[0]) & (fn(xy) < XYRi[1])
    else:
        ind=np.ones_like(xy, dtype=float)

    xy_sub=xy[ind]

    xyP=xyp[0]+1j*xyp[1]
    dxyP=np.diff(xyP)
    lxyP=np.abs(dxyP)
    lP = np.zeros_like(xyP, dtype=float)
    lP[1:]=np.cumsum(lxyP)

    p2 = np.zeros_like(xy_sub, dtype=float)+np.Inf
    atc = np.zeros_like(xy_sub)+np.nan+1j*np.nan

    for xyp_i, dxyp_i, dlp_i, lp_i in zip(xyP[:-1], dxyP, lxyP, lP[:-1]):
        seg_dist_cplx(xy_sub, p2, atc, xyp_i, dxyp_i, lp_i, dlp_i, d_max)

    at_full = np.zeros_like(xy, dtype=float)+np.nan
    at_full[ind]=np.real(atc)
    xt_full = np.zeros_like(xy, dtype=float)+np.nan
    xt_full[ind] = np.imag(atc)
    dist = np.zeros_like(xy, dtype=float)+np.nan
    dist[ind] = np.sqrt(p2)
    return at_full, xt_full, dist
