from numpy import *
from masternodes import masternodes

def localbasis(porder,dim,elemtype):

    pelem,telem,pface,tface,perm = masternodes(porder,dim,elemtype)[0:5]

    npe = pelem.shape[0];
    npf = pface.shape[0];

    if dim==1:
        nve = 2;
        nvf = 1;
    elif dim==2 and elemtype==0: # tri
        nve = 3;
        nvf = 2;
    elif dim==2 and elemtype==1: # quad
        nve = 4;
        nvf = 2;
    elif dim==3 and elemtype==0: # tet
        nve = 4;
        nvf = 3;
    elif dim==3 and elemtype==1: # hex
        nve = 8;
        nvf = 4;

    phielem = zeros((npe, nve));
    phiface = zeros((npf, nvf));

    if dim==1:
        xi  = pelem[:,1-1];
        phielem[:,1-1] = 1.0 - xi;
        phielem[:,2-1] = xi;
        phiface = [1.0];
    elif dim==2 and elemtype==0: # tri
        xi  = pelem[:,1-1];
        eta = pelem[:,2-1];
        phielem[:,1-1] = 1.0 - xi - eta;
        phielem[:,2-1] = xi;
        phielem[:,3-1] = eta;

        xi = pface[:,1-1];
        phiface[:,1-1] = 1.0 - xi;
        phiface[:,2-1] = xi;
    elif dim==2 and elemtype==1: # quad
        xi  = pelem[:,1-1];
        eta = pelem[:,2-1];
        phielem[:,1-1] = (1.0 - xi)*(1.0 - eta);
        phielem[:,2-1] = xi*(1.0 - eta);
        phielem[:,3-1] = xi*eta;
        phielem[:,4-1] = (1.0 - xi)*eta;

        xi = pface[:,1-1];
        phiface[:,1-1] = 1.0 - xi;
        phiface[:,2-1] = xi;
    elif dim==3 and elemtype==0: # tet
        xi   = pelem[:,1-1];
        eta  = pelem[:,2-1];
        zeta = pelem[:,3-1];
        phielem[:,1-1] = 1.0 - xi - eta - zeta;
        phielem[:,2-1] = xi;
        phielem[:,3-1] = eta;
        phielem[:,4-1] = zeta;

        xi = pface[:,1-1];
        eta = pface[:,2-1];
        phiface[:,1-1] = 1.0 - xi - eta;
        phiface[:,2-1] = xi;
        phiface[:,3-1] = eta;
    elif dim==3 and elemtype==1: # hex
        xi   = pelem[:,1-1];
        eta  = pelem[:,2-1];
        zeta = pelem[:,3-1];
        phielem[:,1-1] = (1.0 -xi)*(1.0 -eta)*(1.0 -zeta);
        phielem[:,2-1] = xi*(1.0 -eta)*(1.0 -zeta);
        phielem[:,3-1] = xi*eta*(1.0 -zeta);
        phielem[:,4-1] = (1.0 -xi)*eta*(1.0 -zeta);
        phielem[:,5-1] = (1.0 -xi)*(1.0 -eta)*(zeta);
        phielem[:,6-1] = xi*(1.0 -eta)*(zeta);
        phielem[:,7-1] = xi*eta*(zeta);
        phielem[:,8-1] = (1.0 -xi)*eta*(zeta);

        xi = pface[:,1-1];
        eta = pface[:,2-1];
        phiface[:,1-1] = (1.0 - xi)*(1.0 - eta);
        phiface[:,2-1] = xi*(1.0 - eta);
        phiface[:,3-1] = xi*eta;
        phiface[:,4-1] = (1.0 - xi)*eta;

    return phielem,phiface,pelem,pface,perm
