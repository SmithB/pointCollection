#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 14:44:54 2024

@author: ben
"""

import numpy as np
from .RDE import RDE

def select_valid_pts(D_pt, D_grid, max_dist, max_dt=None):
    '''
    Choose the point data that will be present for any allowable shift

    Parameters
    ----------
    D_pt : pointCollection.data
        Point altimetry data to which the DEM is registered.
    D_grid : pointCollection.grid.data
        grid data to be registered.
    max_dist : float
        maximum distance that the DEM might be displaced.
    max_dt : float
        maximum allowed time between a data point and the DEM

    Returns
    -------
    good : np.array (boolean)
        Points that should be used.

    '''
    if hasattr(D_grid, 't') and D_grid.t is not None:
        good = np.abs(D_pt.t - D_grid.t) < max_dt

    deltas = np.meshgrid(*[np.array([-1, 0, 1])*max_dist for ii in [0, 1]])
    for dx, dy in zip(deltas[0].ravel(), deltas[1].ravel()):
        good[good] &= np.isfinite(D_grid.interp(D_pt.x[good]+dx, D_pt.y[good]+dy, band=0))
    return good

def eval_grid_shift(delta, D_pt, D_grid, sigma_min=0, iterations=1, mask=None):
    '''
    Calculate the misfit for a shift value

    Parameters
    ----------
    delta : iterable (2)
        Shift in x, y.
    D_pt : pointCollection.data
        Reference point data.
    D_grid : pointCollection.grid.data
        grid data to be registered
    sigma_min : float, optional
        Minimum value for data errors. The default is 0.
    iterations : int, optional
        Number of iterations used to discard outliers. The default is 1.
    mask : numpy.array, optional
        Boolean array indicating which values in D_pt to use. The default is None.

    Returns
    -------
    sigma : float
        standard deviation of residuals.
    sigma_scaled : float
        standard deviation of error-scaled residuals.
    mask : numpy.array
        Boolean array indicating which values in D_pt were used.
    m : numpy.array
        Vector indicating mean vertical shift, x-gradient of shift, y_gradient of shift, (if available) rate of shift change
    r : numpy.array
        Residuals to best-fitting model.
    dh : numpy.array
        Differences between point data and DEM.

    '''

    dh = D_pt.z - D_grid.interp(D_pt.x+delta[0], D_pt.y+delta[1], band=0)

    if hasattr(D_grid,'t') and hasattr(D_pt,'t'):
        G=np.c_[np.ones_like(D_pt.x).ravel(),\
                (D_pt.x.ravel() - np.nanmean(D_pt.x)),\
                (D_pt.y.ravel() - np.nanmean(D_pt.y)), \
                (D_pt.t.ravel() - D_grid.t)]
    else:
        G=np.c_[np.ones_like(D_pt.x).ravel(),\
            (D_pt.x.ravel()-np.nanmean(D_pt.x)),\
            (D_pt.y.ravel() - np.nanmean(D_pt.y)) ]

    # Normalize the columns of G to avoid floating-point problems:
    column_scaling = np.ones(G.shape[1])
    for col in range(1, G.shape[1]):
        column_scaling[ col] = np.std(G[:, col])
        G[:, col] /= column_scaling[ col]

    if mask is None:
        mask = np.ones_like(D_pt.x, dtype=bool)

    mask &= np.isfinite(dh)

    last_mask = mask
    iteration=0
    while iteration==0 or (iteration < iterations and not np.all(last_mask==mask)):
        Cinv=(1/(np.maximum(sigma_min, D_pt.sigma[mask]))**2)
        try:
            m=np.linalg.solve(G[mask,:].T @ (np.tile(Cinv[:, None], [1, G.shape[1]])*G[mask,:]), G[mask,:].T @ (Cinv*dh[mask]))
        except np.linalg.LinAlgError:
            m = np.zeros(G.shape[1])
            m[0] = np.nanmean(dh[mask])
            print(f"register_grid: linalg error at iteration {iteration}")
        r=dh-G.dot(m)
        rs = r/np.maximum(sigma_min, D_pt.sigma)
        sigma_hat = RDE(rs[mask])
        last_mask=mask
        mask &= (np.abs(rs) < 3*sigma_hat)
        iteration += 1
    mask=last_mask
    sigma_scaled = np.std(rs[mask])
    sigma = np.std(r[mask])
    return sigma, sigma_scaled, mask, m*column_scaling, r, dh



def search_offsets(D_pt, D_grid, max_delta=10, delta_tol = 1, sigma_min=0.02):
    '''
    Search offsets to find the optimum residual

    Parameters
    ----------
    D_pt : pointCollection.data
        Reference point data.
    D_grid : pointCollection.grid.data
        grid data to be registered.
    max_delta : float, optional
        largest offset attempted. The default is 10.
    delta_tol : float, optional
        Convergence tolerance for offset search. The default is 1.
    sigma_min : float, optional
        Minimum value for estimated data errors. The default is 0.1.

    Returns
    -------
    best_offset : list (2)
        x and y shift for minimum residual.
    R_of_delta : dict
        Dict specifying the residual for all tested offsets.

    '''

    R_of_delta={}
    delta=max_delta/2
    best_offset = [0, 0]
    count=0
    while delta >= delta_tol:
        #print(delta)
        deltas = [ii.ravel() for ii in np.meshgrid(*[np.array([-1, 0, 1])*delta for jj in [0, 1]])]
        for dx, dy  in zip(deltas[0], deltas[1]):
            this_offset = (best_offset[0] + dx, best_offset[1]+dy)
            #print('\t'+str(this_offset))
            if np.any(np.abs(this_offset)> max_delta):
                continue
            if this_offset not in R_of_delta:
                R_of_delta[this_offset] = eval_grid_shift(this_offset, D_pt, D_grid, sigma_min=sigma_min, iterations=1)[1]
                #print('\t'+str([this_offset, R_of_delta[this_offset]]))
        searched=list(R_of_delta.keys())
        Rvals = [R_of_delta[ii] for ii in searched]
        best_offset = searched[np.argmin(Rvals)]
        #print('best_offset:'+str(best_offset))
        shrink = True
        for dx, dy in zip(deltas[0], deltas[1]):
            test_offset = (best_offset[0]+dx, best_offset[1]+dy)
            #print([test_offset, test_offset in R_of_delta])
            if test_offset not in R_of_delta and np.all(np.abs(test_offset) <= max_delta):
                # This case will occur if the best offset is on the boundary of the deltas
                # tested at the current scale, which implies that the search area needs
                # to be shifted rather than shrunk
                shrink=False
                break
        if shrink:
            delta /= 2
        count += 1
        if count > 20:
            print('quitting because count is too large\n')
            break

    return best_offset, R_of_delta
