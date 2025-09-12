from numpy import *
#from koornwinder import *

import importlib
import koornwinder
importlib.reload(koornwinder)
from koornwinder import *

def tensorproduct(x,porder):

    if len(x.shape)==1:
        n = x.shape[0]
        x = reshape(x,(n,1))
        dim = 1
    else:
        dim = x.shape[1]
        n = x.shape[0]

    if dim==1:
        f,fx,fy,fz = koornwinder(x,porder)
    elif dim==2:
        g1,gx = koornwinder(x[:,0],porder)[0:2]
        g2,gy = koornwinder(x[:,1],porder)[0:2]
        f  = zeros((n,(porder+1)*(porder+1)))
        fx = zeros((n,(porder+1)*(porder+1)))
        fy = zeros((n,(porder+1)*(porder+1)))
        fz = 0
        for ii in range(0,n):
            f[ii,:] =  kron(g2[ii,:],g1[ii,:])
            fx[ii,:] = kron(g2[ii,:],gx[ii,:])
            fy[ii,:] = kron(gy[ii,:],g1[ii,:])
    elif dim==3:
        g1,gx = koornwinder(x[:,0],porder)[0:2]
        g2,gy = koornwinder(x[:,1],porder)[0:2]
        g3,gz = koornwinder(x[:,2],porder)[0:2]
        f  = zeros((n,(porder+1)*(porder+1)*(porder+1)))
        fx = zeros((n,(porder+1)*(porder+1)*(porder+1)))
        fy = zeros((n,(porder+1)*(porder+1)*(porder+1)))
        fz = zeros((n,(porder+1)*(porder+1)*(porder+1)))
        for ii in range(0,n):
            f[ii,:] =  kron(g3[ii,:],kron(g2[ii,:],g1[ii,:]))
            fx[ii,:] = kron(g3[ii,:],kron(g2[ii,:],gx[ii,:]))
            fy[ii,:] = kron(g3[ii,:],kron(gy[ii,:],g1[ii,:]))
            fz[ii,:] = kron(gz[ii,:],kron(g2[ii,:],g1[ii,:]))

    return f,fx,fy,fz
