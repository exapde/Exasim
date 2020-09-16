from numpy import *
from mkelconcg import mkelconcg
from mkent2elem import mkent2elem
from xiny import xiny

def mkcgent2dgent(dgnodes,tol):

    dig = log10(1/tol);
    dgnodes = dgnodes.round(int(dig));

    # CG element-to-node connectivities
    cgnodes,cgelcon = mkelconcg(dgnodes)[0:2];

    # CG node-to-element connectivities
    rowent2elem,colent2elem = mkent2elem(cgelcon)[0:2];

    cgent2dgent = colent2elem.flatten('F');
    npe = dgnodes.shape[0];
    dim = dgnodes.shape[1];
    nent = max(cgelcon.flatten())+1;
    for i in range(0,nent): # for CG node i
        nelem = rowent2elem[i+1]-rowent2elem[i]; # number of elements connected to CG node i
        elem = colent2elem[rowent2elem[i]:rowent2elem[i+1]]; # elements connected to CG node i
        xcg = reshape(cgnodes[:,i],(1,dim),'F'); # coordinates of CG node i
        for j in range(0,nelem): # for element j connected to CG node i
            xdg = dgnodes[:,:,elem[j]]; # coordinates of DG nodes on element j
            ind = xiny(xcg,xdg,1); # match CG node i to one of the DG nodes on element j
            cgent2dgent[rowent2elem[i]+j] = elem[j]*npe+ind[0]; # index of the matched DG node

    return cgelcon,rowent2elem,colent2elem,cgent2dgent,cgnodes
