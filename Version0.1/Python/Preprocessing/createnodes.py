from numpy import *
from masternodes import masternodes

def createnodes(p,t,porder):

    nve = t.shape[0];
    ne = t.shape[1];
    nd = p.shape[0];

    elemtype = 0;
    if ((nd==2) and (nve==4)):
        elemtype=1;
    if ((nd==3) and (nve==8)):
        elemtype=1;

#    print([porder,nd,elemtype])
    plocal,tm1,tm2,tm3,perm = masternodes(porder,nd,elemtype)[0:5];

    npl = plocal.shape[0]
    nt = ne
    ns = nve

    philocal = zeros((npl,ns))
    if nd==1:
        xi  = plocal[:,0]
        philocal[:,0] = 1.0 - xi
        philocal[:,1] = xi
    elif nd==2 and ns==3:
        xi  = plocal[:,0]
        eta = plocal[:,1]
        philocal[:,0] = 1.0 - xi - eta
        philocal[:,1] = xi
        philocal[:,2] = eta
    elif nd==2 and ns==4:
        xi  = plocal[:,0]
        eta = plocal[:,1]
        philocal[:,0] = (1-xi)*(1-eta)
        philocal[:,1] = xi*(1-eta)
        philocal[:,2] = xi*eta
        philocal[:,3] = (1-xi)*eta
    elif nd==3 and ns==4:
        xi   = plocal[:,0]
        eta  = plocal[:,1]
        zeta = plocal[:,2]
        philocal[:,0] = 1 - xi - eta - zeta
        philocal[:,1] = xi
        philocal[:,2] = eta
        philocal[:,3] = zeta
    elif nd==3 and ns==8:
        xi   = plocal[:,0]
        eta  = plocal[:,1]
        zeta = plocal[:,2]
        philocal[:,0] = (1-xi)*(1-eta)*(1-zeta)
        philocal[:,1] = xi*(1-eta)*(1-zeta)
        philocal[:,2] = xi*eta*(1-zeta)
        philocal[:,3] = (1-xi)*eta*(1-zeta)
        philocal[:,4] = (1-xi)*(1-eta)*(zeta)
        philocal[:,5] = xi*(1-eta)*(zeta)
        philocal[:,6] = xi*eta*(zeta)
        philocal[:,7] = (1-xi)*eta*(zeta)

    dgnodes=zeros((npl,nt,nd))
    for dim in range(0,nd):
        for node in range(0,ns):
            tm1 = reshape(philocal[:,node],(npl,1))
            tm2 = p[dim,t[node,:]]
            tm2 = reshape(tm2,(1,len(tm2)))
            dp = dot(tm1,tm2)
            dgnodes[:,:,dim] = dgnodes[:,:,dim] + dp

    dgnodes = dgnodes.swapaxes(1,2)

    return dgnodes, elemtype, perm
