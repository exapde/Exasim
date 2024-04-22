from numpy import *
# from tensorproduct import *
# from koornwinder import *

import importlib
import tensorproduct
import koornwinder
importlib.reload(tensorproduct)
importlib.reload(koornwinder)
from tensorproduct import *
from koornwinder import *

def mkshape(porder,plocal,pts,elemtype):

    if porder==0:
        d = plocal.size
        m = pts.shape[0]
        nfs = zeros((1,m,d+1))
        nfs[:,:,0] = 1.0
    else:
        if elemtype==0:
            nf,nfx,nfy,nfz = koornwinder(pts,porder)
            A = koornwinder(plocal,porder)[0]
        else:
            nf,nfx,nfy,nfz = tensorproduct(pts,porder)
            A = tensorproduct(plocal,porder)[0]
        Ainv =  linalg.inv(A)

        m = nf.shape[0]
        n = A.shape[0]
        d = plocal.shape[1]
        nfs = zeros((n,m,d+1))
        if d==1:
            tm = dot(nf,Ainv)
            nfs[:,:,0] = tm.transpose()
            tm = dot(nfx,Ainv)
            nfs[:,:,1] = tm.transpose()
        elif d==2:
            tm = dot(nf,Ainv)
            nfs[:,:,0] = tm.transpose()
            tm = dot(nfx,Ainv)
            nfs[:,:,1] = tm.transpose()
            tm = dot(nfy,Ainv)
            nfs[:,:,2] = tm.transpose()
        elif d==3:
            tm = dot(nf,Ainv)
            nfs[:,:,0] = tm.transpose()
            tm = dot(nfx,Ainv)
            nfs[:,:,1] = tm.transpose()
            tm = dot(nfy,Ainv)
            nfs[:,:,2] = tm.transpose()
            tm = dot(nfz,Ainv)
            nfs[:,:,3] = tm.transpose()

    return nfs
