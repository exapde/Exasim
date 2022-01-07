from numpy import *
from localbasis import localbasis
from getelemface import getelemface
from xiny import xiny
import sys

def faceconnectivity(p,t,f,f2t,elemtype,prdexpr,porder):

    dim = p.shape[0];
    nprd = len(prdexpr);
    nf = f2t.shape[1];

    philocvl,philocfc,tm1,tm2,perm = localbasis(porder,dim,elemtype)[0:5]
    philocvl = philocvl.transpose();
    philocfc = philocfc.transpose();
    perm = perm-1;
    face = getelemface(dim,elemtype);

    npf = perm.shape[0];
    npe = philocvl.shape[1];

    facecon = zeros((npf,2,nf)).astype(int);
    for i in range(0,nf):
        e1 = f2t[0,i]-1;
        l1 = f2t[1,i]-1;
        e2 = f2t[2,i]-1;
        l2 = f2t[3,i]-1;

        f1 = t[face[:,l1],e1];
        pf = dot(p[:,f1.flatten('F')],philocfc); # nodal points on face i
        pf = pf.round(8);

        t1 = t[:,e1];
        p1 = dot(p[:,t1],philocvl);
        p1 = p1[:,perm[:,l1].flatten('F')];
        p1 = p1.round(8);

        if (f[l1,e1]<0): # periodic face
            pbnd = -f[l1,e1];
            for k in range(0,nprd):
                if (prdexpr[k,0]==pbnd):
                    pf = prdexpr[k,1](pf);
                    p1 = prdexpr[k,1](p1);
                    break;
                elif (prdexpr[k,2]==pbnd):
                    pf = prdexpr[k,3](pf);
                    p1 = prdexpr[k,3](p1);
                    break;

        j1 = xiny(p1,pf,0);
        facecon[:,0,i] = e1*npe + perm[j1,l1];

        if (e2>=0): # face i is an interior face
            t2 = t[:,e2];
            p2 = dot(p[:,t2],philocvl);
            p2 = p2[:,perm[:,l2].flatten('F')];
            p2 = p2.round(8);
            if (f[l2,e2]<0): # periodic face
                pbnd = -f[l2,e2];
                if (prdexpr[k,0]==pbnd):
                    p2 = prdexpr[k,1](p2);
                elif (prdexpr[k,2]==pbnd):
                    p2 = prdexpr[k,3](p2);
            j2 = xiny(p2,pf,0);
            facecon[:,1,i] = e2*npe + perm[j2,l2];
        else:
            facecon[:,1,i] = facecon[:,0,i];

    #print(reshape(facecon[:,0,:],(npf,nf),'F').T)

    return facecon
