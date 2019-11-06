# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 15:44:31 2018

@author: ben
"""

import numpy as np
def pt_blockmedian(xx, yy, zz, delta, xy0=[0.,0.], return_index=False):
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
    xm=np.zeros_like(ux)+np.NaN
    ym=np.zeros_like(ux)+np.NaN
    zm=np.zeros_like(ux)+np.NaN
    if return_index:
        ind=np.zeros((ux.size,2), dtype=int)
    
    ix=np.hstack((ix.ravel(), xyind.size))
    for  count, i0 in enumerate(ix[:-1]):
        ii=np.arange(i0, ix[count+1], dtype=int)
        iM=np.maximum(ii.size/2.-1, 0)
        if (iM-np.floor(iM)==0) and ii.size>1:
            iM=int(np.floor(iM))
            iM=ii[[iM, iM+1]]
            xm[count]=(xs[iM[0]]+xs[iM[1]])/2.
            ym[count]=(ys[iM[0]]+ys[iM[1]])/2.
            zm[count]=(zs[iM[0]]+zs[iM[1]])/2.
            if return_index:
                ind[count,:]=iM
        else:
            try:
                iM=ii[int(iM)]
            except IndexError:
                print("HERE")
            xm[count]=xs[iM]
            ym[count]=ys[iM]
            zm[count]=zs[iM]
            if return_index:
                ind[count,:]=iM
        #plt.figure(1); plt.clf(); plt.subplot(211); plt.cla(); plt.scatter(xs[ii]-xm[count], ys[ii]-ym[count], c=zs[ii]); plt.axis('equal'); plt.subplot(212); plt.plot(zs[ii]); 
        #plt.pause(1)
        #print(count)
    if return_index:
        return xm, ym, zm, ind
    else:
        return xm, ym, zm

