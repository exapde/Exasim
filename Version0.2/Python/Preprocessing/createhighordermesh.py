import numpy as np
from createnodes import createnodes
from facenumbering import facenumbering

def createhighordermesh(mesh,app):

    print("run createnodes...");
    mesh['dgnodes'], elemtype, perm = createnodes(mesh['p'],mesh['t'],app['porder'])[0:3];

    print("run facenumbering...");
    mesh['f'], mesh['tprd'] = facenumbering(mesh['p'],mesh['t'],elemtype,mesh['boundaryexpr'],mesh['periodicexpr'])[0:2];

    # Project nodes on the curved boundaries
    if (len(mesh['curvedboundaryexpr'])>0) and (app['porder']>1) and (max(mesh['curvedboundary'])>0):
        print("Project dgnodes onto the curved boundaries...");
        fd = mesh['curvedboundaryexpr'];
        nd = shape(mesh.p,0);
        nfe,ne = shape(mesh.f);
        if (nd == 2):
            for i in range(0,ne):
                for j in range(0,nfe):
                    if (mesh['f'][j,i]!=0): # boundary element
                        k = abs(mesh.f[j,i]); # get boundary index
                        if (mesh['curvedboundary'][k]==1): # if this boundary is curved
                            p = mesh['dgnodes'][perm[:,j],:,i];
                            deps = 1e-8*max(max(p,axis=0)-min(p,axis=0));
                            d = fd[k](p);
                            dgradx = (fd[k]([p[:,0]+deps, p[:,1]])-d)/deps;
                            dgrady = (fd[k]([p[:,0], p[:,1]+deps])-d)/deps;
                            dgrad2 = dgradx*dgradx + dgrady*dgrady;
                            dgrad2[dgrad2==0] = 1;
                            p = p-[d*dgradx/dgrad2, d*dgrady/dgrad2];
                            mesh['dgnodes'][perm[:,j],:,i] = p;
        elif (nd==3):
            for i in range(0,ne):
                for j in range(0,nfe):
                    if (mesh['f'][j,i]!=0): # boundary element
                        k = abs(mesh.f[j,i]); # get boundary index
                        if (mesh['curvedboundary'][k]==1): # if this boundary is curved
                            p = mesh['dgnodes'][perm[:,j],:,i];
                            deps = 1e-8*max(max(p,axis=0)-min(p,axis=0));
                            d = fd[k](p);
                            dgradx = (fd[k]([p[:,0]+deps,p[:,1],p[:,2]])-d)/deps;
                            dgrady = (fd[k]([p[:,0],p[:,1]+deps,p[:,2]])-d)/deps;
                            dgradz = (fd[k]([p[:,0],p[:,1],p[:,2]+deps])-d)/deps;
                            dgrad2 = dgradx*dgradx+dgrady*dgrady+dgradz*dgradz;
                            dgrad2[dgrad2==0] = 1;
                            p = p-[d*dgradx/dgrad2, d*dgrady/dgrad2, d*dgradz/dgrad2];
                            mesh['dgnodes'][perm[:,j],:,i] = p;

    return mesh
