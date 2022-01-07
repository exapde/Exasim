from numpy import *

def  createcgcells(cgelcon,tlocal,ne):

    nce, nve = tlocal.shape;
    cells = zeros((ne*nce,nve));
    tlocal = tlocal.flatten('F')-1;
    for el in range (0,ne):
        m = nce*el;
        cells[m:(m+nce),:] = reshape(cgelcon[tlocal,el],[nce,nve],'F');
    #cells = cells - 1;

    return cells
