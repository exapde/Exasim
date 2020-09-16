from numpy import *
from mkcgent2dgent import mkcgent2dgent
from createcgcells import createcgcells
from getcelltype import getcelltype

def createcggrid(dgnodes,tlocal):

    # Get some dimensions of the mesh
    nd = dgnodes.shape[1];
    ne = dgnodes.shape[2];
    nve = tlocal.shape[1];

    cgelcon,tm1,tm2,tm3,cgnodes = mkcgent2dgent(dgnodes,1e-8)[0:5];

    cgnodes = cgnodes.T;
    cgcells = createcgcells(cgelcon,tlocal,ne);

    elemtype = 0;
    if (nd==2) and (nve==4):
        elemtype=1;
    if (nd==3) and (nve==8):
        elemtype=1;

    cell_t = getcelltype(nd,elemtype);

    return cgnodes, cgelcon, cgcells, cell_t
