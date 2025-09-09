import numpy as np
import scipy
import pointCollection as pc
def calc_bias_and_R(D, i0, i1, R, min_sigma=0.5):
    x0=D.x[i0]
    y0=D.y[i0]
    z0=D.z[i0]
    i1 = i1[np.flatnonzero(~np.in1d(i1, i0))]
    dx=D.x[i1]-x0
    dy=D.y[i1]-y0
    dz=D.z[i1]-z0
    ii=np.ones_like(i1, dtype=bool)
    
    constraint_eqs=0.001*np.array([[1, 0, 0], [0, 1, 0]])
    constraints=np.array([0., 0.])
    G=np.c_[dx/1000, dy/1000, np.ones_like(dx)]
    
    count=0
    while np.sum(ii)> 4:
        count +=1
        G1=np.vstack((G[ii,:], constraint_eqs))
        d1=np.concatenate([dz[ii], constraints], axis=0)
        m=scipy.linalg.lstsq(G1, d1, lapack_driver='gelsy')[0]
        r=dz-G.dot(m)
        sigma_hat=pc.RDE(r[ii])
        if count > 5:
            break
        last_ii=ii.copy()
        ii = np.abs(r) < 3*np.maximum(sigma_hat, min_sigma)
        if np.all(last_ii == ii):
            break
    
    N=np.sum(ii)
    # the negative sign here is legacy of the Matlab version of the code.
    bias=-m[-1]
    noise=np.std(r[ii])
    slope=[m[0]*1000, m[1]*1000]
    
    return bias, noise, slope, N
