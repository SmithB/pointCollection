# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 15:44:31 2018

@author: ben
"""

import numpy as np
def pt_blockmedian(xx, yy, zz, delta, xy0=[0.,0.], return_index=False, index_only=False, index_and_count_only=False, return_count=False):

    if index_only or index_and_count_only or return_count:
       return_index=True
    if index_and_count_only:
        return_count=True
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
    sorted_ind=np.argsort(xyind+(z-np.nanmin(z))/zscale)
    xs=x[sorted_ind]
    ys=y[sorted_ind]
    zs=z[sorted_ind]
    xyind=xyind[sorted_ind]
    ux, ix=np.unique(xyind, return_index=True)
    if not index_only:
        xm=np.zeros_like(ux)+np.NaN
        ym=np.zeros_like(ux)+np.NaN
        zm=np.zeros_like(ux)+np.NaN
    if return_index:
        ind=np.zeros((ux.size,2), dtype=int)
    if return_count:
        N=np.zeros(ux.size, dtype=int)
        
    ix=np.hstack((ix.ravel(), xyind.size))
    for  count, i0 in enumerate(ix[:-1]):
        # number of elements in bin
        ni=ix[count+1]-i0
        if return_count:
            N[count]=ni
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
                ind[count,:]=sorted_ind[iM]
        else:
            # odd case: iM ends in 0.5.  if ni=3, iM=3/2-1 = 0.5, want 1
            # if ni=1, iM=-0.5, want 0
            iM=int(i0+np.floor(iM)+1)
            if not index_only:
                xm[count]=xs[iM]
                ym[count]=ys[iM]
                zm[count]=zs[iM]
            if return_index:
                ind[count,:]=sorted_ind[iM]
        #plt.figure(1); plt.clf(); plt.subplot(211); plt.cla(); plt.scatter(xs[ii]-xm[count], ys[ii]-ym[count], c=zs[ii]); plt.axis('equal'); plt.subplot(212); plt.plot(zs[ii]); 
        #plt.pause(1)
        #print(count)
    if index_only:
        return ind
    if index_and_count_only:
        return ind, N
    if return_index:
        if return_count:
            return xm, ym, zm, ind, N
        else:
            return xm, ym, zm, ind 
    else:
        if return_count:
            return xm, ym, zm, N
        else:
            return xm, ym, zm,


def __main__():
    scale=10
    N=200000
    delta=2
    x=np.random.rand(N)*scale
    y=np.random.rand(N)*scale
    z=np.random.rand(N)*scale

    xm, ym, zm0, ii=pt_blockmedian(x, y, z, delta, return_index=True)

    zm=0.5*(z[ii[:,0]]+z[ii[:,1]])

    xi=np.floor(x/delta)*delta
    yi=np.floor(y/delta)*delta
    if False:
        for x0 in np.unique(xi):
            for y0 in np.unique(yi):
                els=(xi==x0) & (yi==y0)
                if ~np.any(els):
                    continue
                zmi=np.median(z[els])

                this=(np.floor(xm/delta)*delta == x0) & (np.floor(ym/delta)*delta == y0)
                print(f"({x0}, {y0}: {zmi-zm[this]}?")

if __name__=='__main__':
    __main__()

