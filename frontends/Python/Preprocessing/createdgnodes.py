from numpy import *
from masternodes import masternodes

def createdgnodes(p,t,f,curvedboundary,curvedboundaryexpr,porder):

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

    # Project nodes on the curved boundaries
    if (len(curvedboundaryexpr)>0) and (porder>1) and (max(curvedboundary)>0):
        print("Project dgnodes onto the curved boundaries...");
        fd = curvedboundaryexpr;
        nd = p.shape[0];
        nfe,ne = f.shape;
        perm = perm-1;
        if (nd == 2):
            for i in range(0,ne):
                for j in range(0,nfe):
                    if (f[j,i]!=0): # boundary element
                        k = abs(f[j,i])-1; # get boundary index
                        if (curvedboundary[k]==1): # if this boundary is curved
                            p = dgnodes[perm[:,j],:,i];
                            deps = 1e-8*amax(amax(p,axis=0)-amin(p,axis=0));
                            d = fd[k](p.T);
                            tm = row_stack((p[:,0]+deps, p[:,1]));
                            dgradx = (fd[k](tm)-d)/deps;
                            tm = row_stack((p[:,0], p[:,1]+deps));
                            dgrady = (fd[k](tm)-d)/deps;
                            dgrad2 = dgradx*dgradx + dgrady*dgrady;
                            dgrad2[dgrad2==0] = 1;
                            p = p - column_stack((d*dgradx/dgrad2, d*dgrady/dgrad2));
                            dgnodes[perm[:,j],:,i] = p;
                            
        elif (nd==3):
            for i in range(0,ne):
                for j in range(0,nfe):
                    if (f[j,i]!=0): # boundary element
                        k = abs(f[j,i])-1; # get boundary index
                        if (curvedboundary[k]==1): # if this boundary is curved
                            p = dgnodes[perm[:,j],:,i];
                            deps = 1e-8*max(max(p,axis=0)-min(p,axis=0));
                            d = fd[k](p.T);
                            tm = row_stack((p[:,0]+deps,p[:,1],p[:,2]))
                            dgradx = (fd[k](tm)-d)/deps;
                            tm = row_stack((p[:,0],p[:,1]+deps,p[:,2]))
                            dgrady = (fd[k](tm)-d)/deps;
                            tm = row_stack((p[:,0],p[:,1],p[:,2]+deps))
                            dgradz = (fd[k](tm)-d)/deps;
                            dgrad2 = dgradx*dgradx+dgrady*dgrady+dgradz*dgradz;
                            dgrad2[dgrad2==0] = 1;
                            p = p - column_stack((d*dgradx/dgrad2, d*dgrady/dgrad2, d*dgradz/dgrad2));
                            dgnodes[perm[:,j],:,i] = p;

    return dgnodes
