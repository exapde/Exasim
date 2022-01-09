from numpy import *
#from importlib import reload

def mkelemblocks(ne,ns):

    ns = min(ns,ne);
    nb = floor(ne/ns);  # number of blocks
    na = round(ne/nb); # number of elements per block
    nk = arange(1,ne,na);
    m = nk.size;
    if (m==1):
        nm = array([1, ne]).reshape(2,1);
    else:
        nm = zeros((2,m))
        nm[0,:] = nk;
        nm[1,:] = concatenate([nk[1:]-1, [ne]]);
        if ((nm[1,-1]-nm[0,-1]) < (na/2)):
            nm[1,-2] = nm[1,-1];
            nm = nm[:,0:-1]

    nb = nm.shape[1];
    while (min(nm[1,:]-nm[0,:])<(ns/2-1) or max(nm[1,:]-nm[0,:])>ns):
        nb = nb+1;
        na = round(ne/nb); # number of elements per block
        nk = arange(1,ne,na);
        m = nk.size;
        nm = zeros((2,m))
        nm[0,:] = nk;
        nm[1,:] = concatenate([nk[1:]-1, [ne]]);
        if (nm[1,-1]-nm[0,-1]) < (na/2):
            nm[1,-2] = nm[1,-1];
            nm = nm[:,0:-1]
    nb = nm.shape[1];

    if nm[1,-1]!=ne:
        error("something wrong");

    if max(nm[1,:]-nm[0,:])>ns:
        error("something wrong");

    if min(nm[1,:]-nm[0,:])<(ns/2-1):
        error("something wrong");

    return nm,nb
