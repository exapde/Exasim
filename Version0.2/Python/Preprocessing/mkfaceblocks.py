from numpy import *
from mkelemblocks import mkelemblocks

def mkfaceblocks(mf,bcm,ns):

    for i in range(0,mf.size-1):
        nf = mf[i+1]-mf[i]
        nmf,nb = mkelemblocks(nf,ns);
        tm = vstack((mf[i]+nmf,bcm[i]*ones((1,nmf.shape[1]))));
        if i == 0:
            nm = tm;
        else:
            nm = hstack((nm, tm));
    nb = nm.shape[1];

    return nm,nb
