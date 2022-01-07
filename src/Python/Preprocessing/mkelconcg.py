from numpy import *
from xiny import xiny

def mkelconcg(dgnodes):

    # remove duplicate nodes in mesh.p1
    ns = dgnodes.shape[0];
    dim = dgnodes.shape[1];
    nt = dgnodes.shape[2];
    A = reshape(dgnodes.transpose(0,2,1),(ns*nt,dim),'F');

    # B = unique(A,axis=0);
    # b = xiny(A,B,1);
    B, ia, b = unique(A,return_index=True,return_inverse=True,axis=0);

    # CG mesh
    cgnodes = B.transpose();
    cgelcon = reshape(b, (ns, nt), order='F');

    return cgnodes,cgelcon
