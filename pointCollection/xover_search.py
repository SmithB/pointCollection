# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 15:42:22 2019

@author: ben
"""
import numpy as np
import pointCollection as pc
#import matplotlib.pyplot as plt

def dilate_bins(bins, delta):
    '''
    to a list of bins, add all neighboring bins
    '''
    original_bins=bins.copy()
    for delta_x in [-1, 0, 1]:
        for delta_y in [-1, 0, 1]:
            for this_bin in original_bins:
                this_shifted_bin=(this_bin[0]+delta_x*delta, this_bin[1]+delta_y*delta)
                bins.add(this_shifted_bin)
 
def x_point(A, B):
    '''
    find the point at which vectors A and B cross
    '''

    dA=A[-1]-A[0]
    dB=B[-1]-B[0]
    det=-np.imag(dA*(dB.conjugate()))
    if det==0:
        return None, None
    dAB0=A[0]-B[0]
    lA=np.imag(dAB0*(dB.conjugate()))/det
    lB=np.imag(dAB0*(dA.conjugate()))/det

    if (lA <0 ) or (lA >1) or (lB<0) or (lB >1):
        return None, None

    return A[0]+lA*dA, [lA, lB]

def cross_by_time(t, xy, ep, l0):
        
    Wsub=[[],[]]
    delta_list=[[0],[0]]
    for ii in [0, 1]:
        t_x=t[ii][ep[ii][0]]+(t[ii][ep[ii][1]]-t[ii][ep[ii][0]])*l0[ii]
        Wsub[ii]=np.where( t[ii] < t_x )[0][-1]
        if (Wsub[ii].size==0) or (Wsub[ii]>=len(t[ii])-1):
            Wsub[ii]=len(t[ii])-2
        if (Wsub[ii] > 1) and ((t_x - t[ii][Wsub[ii]]) < 0.1*(t[ii][Wsub[ii]] - t[ii][Wsub[ii]-1])):
            delta_list[ii] += [-1]
        if (Wsub[ii] < t[ii].size-3) and ((t_x - t[ii][Wsub[ii]]) > 0.9*(t[ii][Wsub[ii]+1] - t[ii][Wsub[ii]])):
            delta_list[ii] += [1]
    for delta12 in [(a,b) for a in delta_list[0] for b in delta_list[1]]:
        try:
            xyC, l0 = x_point(xy[0][[ Wsub[0]+delta12[0], Wsub[0]+1+delta12[0]]], xy[1][[ Wsub[1]+delta12[1], Wsub[1]+1+delta12[1] ]] )
        except IndexError:
            print("Here")
        if xyC is not None:
            Wsub[0] += delta12[0]
            Wsub[1] += delta12[1]
            break
    if xyC is not None:
        ep = [[X, X+1] for X in Wsub]
    return xyC, ep, l0

def cross_by_zoom(T, inds, delta, DOPLOT=False):
    '''
    find the point at which the paths in T cross
    '''
    if inds is None:
        inds=[np.arange(len(Ti.x)) for Ti in T]
    
    xy=[T[ii].x[inds[ii]]+1j*T[ii].y[inds[ii]] for ii in range(len(T))]
    use_time=False
    if hasattr(T[0], 'time'):
        t=[T[ii].time[inds[ii]] for ii in range(len(T))]
        use_time=True
    try:
        xyC, l0=x_point(xy[0][[0, -1]], xy[1][[0, -1]])
    except ValueError:
        print("HERE")
    if xyC is None:
        return None, None, None
    Lsearch=np.nanmax([np.nanmax(np.abs(xyi-xyC)) for xyi in xy])
    xy_ep=[[], []]
    ep=[[ii[0], ii[-1]] for ii  in inds]
    mask=[np.ones_like(xyi, dtype=bool) for xyi in xy]
    while np.max([np.diff(ii) for ii in ep]) > 1:
        ep_last=ep.copy()
        for ii in [0, 1]:
            mask_temp=mask[ii] & (np.abs(xy[ii]-xyC) <= Lsearch)
            if np.sum(mask_temp)<2:
                return None, None, None
            if np.sum(mask_temp)>=2:
                mask[ii]=mask_temp
                W=np.where(mask[ii])[0]
                ep[ii]=[W[0], W[-1]]
                xy_ep[ii]=xy[ii][ep[ii]]    
        if use_time:
            xyC, ep, l0 = cross_by_time(t, xy, ep, l0)
        if ep==ep_last and (Lsearch==delta):
            return None, None, None
        xyC, l0 = x_point(xy[0][[ep[0][0], ep[0][1]]], xy[1][[ep[1][0], ep[1][1]]] )
        if DOPLOT:
            import matplotlib.pyplot as plt
            plt.plot(xyC.real, xyC.imag,'ko')
        if xyC is None:
            return None, None, None
        Lsearch=np.maximum(Lsearch/2, delta)
    return [xyC.real, xyC.imag], [inds[0][ep[0]], inds[1][ep[1]]], l0
        
 
def cross_tracks(T, delta_coarse=10, delta=1, DOPLOT=False):
    '''
    find the point at which the paths in T cross
    '''     
    
    # find the bins that are in both elements of T
    B_fun=[np.round(np.c_[Ti.x, Ti.y]/delta_coarse)*delta_coarse for Ti in T]
    B=[pc.unique_by_rows(Bsub, return_dict=True)[1] for Bsub in B_fun]
    bin_keys=B[0].keys() & B[1].keys()
    if len(bin_keys)==0:
        return None, None, None
    dilate_bins(bin_keys, delta_coarse)
    intersect_inds=[np.array(sorted(np.concatenate([Bi[kk] for kk in bin_keys if kk in Bi]))) for Bi in B]
    if DOPLOT:
        import matplotlib.pyplot as plt
        for ii in [0, 1]:
            plt.plot(T[ii].x[intersect_inds[ii]], T[ii].y[intersect_inds[ii]],'.')
        for BB in bin_keys:
            plt.plot(BB[0], BB[1],'k*')
    xyC, inds, L = cross_by_zoom(T, intersect_inds, delta)
    return xyC, inds, L

def resample_path(x, y, spacing):
    """
    interpolate a path to a different resolution
    
    inputs:
        x, y : path coordinates
        spacing: intended resolution of the result
        
    Calculates the along-track distance for each point in x and y, then interpolates
    the x and y coordinates to an evenly spaced vector of distance values, separated by
    'spacing'
    """     
    
    s0=np.concatenate([[0], np.cumsum(np.diff(np.abs(x+1j*y)))])
    s=np.arange(s0[0], s0[-1]+spacing, spacing)
    if s[-1] < s0[-1]:
        s=np.concatenate([s, [s0[-1]]])
    return np.interp(s, s0, x), np.interp(s, s0, y)


def __test__():
    import matplotlib.pyplot as plt
    import pointCollection as pc
    x0=np.arange(0, 13, .2)
    y0=0.01*(x0*2)**2
    x0, y0=resample_path(x0, y0, 0.1)
    x1=np.arange(0.5, 4.95, .1)
    y1=-0.25*(x1**2)+x1+2
    x1, y1=resample_path(x1, y1, 0.1)
    
    
    plt.figure()
    plt.plot(x0, y0)
    plt.plot(x1, y1)
    plt.axis('equal')
    D=[pc.data().from_dict({'x':x0, 'y':y0, 'time':np.arange(len(x0))}),
       pc.data().from_dict({'x':x1, 'y':y1, 'time':np.arange(len(x1))})]


    xyC, inds, L=cross_tracks(D, delta=0.1, delta_coarse=0.5)

