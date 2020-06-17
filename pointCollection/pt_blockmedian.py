# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 15:44:31 2018

@author: ben
"""

import numpy as np
def pt_blockmedian(xx, yy, zz, delta, xy0=[0.,0.], return_index=False, index_only=False):

    if index_only:
        return_index=True
    temp=np.isfinite(zz)
    x=xx[temp]
    y=yy[temp]
    z=zz[temp]
    if len(x)==0:
        if return_index:
            return np.zeros([0]), np.zeros([0]), np.zeros([0]), np.zeros([0])
        else:
            return np.zeros([0]), np.zeros([0]), np.zeros([0])
    yscale=np.ceil((np.nanmax(y)-np.nanmin(y))/delta*1.1)
    zscale=(np.nanmax(z)-np.nanmin(z))*1.1
    xr=np.floor((x-xy0[0])/delta)
    yr=np.floor((y-xy0[1])/delta)
    xyind=xr*yscale+(yr-np.min(yr))
    ind=np.argsort(xyind+(z-np.nanmin(z))/zscale)
    xs=x[ind]
    ys=y[ind]
    zs=z[ind]
    xyind=xyind[ind]
    ux, ix=np.unique(xyind, return_index=True)
    if not index_only:
        xm=np.zeros_like(ux)+np.NaN
        ym=np.zeros_like(ux)+np.NaN
        zm=np.zeros_like(ux)+np.NaN
    if return_index:
        ind=np.zeros((ux.size,2), dtype=int)
    
    ix=np.hstack((ix.ravel(), xyind.size))
    for  count, i0 in enumerate(ix[:-1]):
        # number of elements in bin
        ni=ix[count+1]-i0
        # index of median in the bin
        iM=ni/2.-1
        if (iM-np.floor(iM)==0) and ni>1:
            iM=int(np.floor(iM))
            iM=np.array([i0+iM, i0+iM+1])
            if not index_only:
                xm[count]=(xs[iM[0]]+xs[iM[1]])/2.
                ym[count]=(ys[iM[0]]+ys[iM[1]])/2.
                zm[count]=(zs[iM[0]]+zs[iM[1]])/2.
            if return_index:
                ind[count,:]=iM
        else:
            # odd case: iM ends in 0.5.  if ni=3, iM=3/2-1 = 0.5, want 1
            # if ni=1, iM=-0.5, want 0
            iM=int(i0+np.floor(iM)+1)
            if not index_only:
                xm[count]=xs[iM]
                ym[count]=ys[iM]
                zm[count]=zs[iM]
            if return_index:
                ind[count,:]=iM
        #plt.figure(1); plt.clf(); plt.subplot(211); plt.cla(); plt.scatter(xs[ii]-xm[count], ys[ii]-ym[count], c=zs[ii]); plt.axis('equal'); plt.subplot(212); plt.plot(zs[ii]); 
        #plt.pause(1)
        #print(count)
    if index_only:
        return ind
    if return_index:
        return xm, ym, zm, ind
    else:
        return xm, ym, zm

