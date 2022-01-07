from numpy import *
from getelemface import getelemface
from xiny import xiny
#from mkf2t import mkf2t
from mkf2e import mkf2e

def  facenumbering(p,t,elemtype,bndexpr,prdexpr):

    dim = p.shape[0]
    nve = t.shape[0]
    ne = t.shape[1]
    face = getelemface(dim,elemtype);
    nvf,nfe = shape(face);

    t2fl = reshape(t[face.flatten(order='F'),:], (nvf, nfe*ne), order='F');
    pf = reshape(p[:,t2fl.flatten(order='F')], (dim, nvf, nfe*ne), order='F');
    pf = reshape(sum(pf,axis=1)/nvf,(dim,nfe,ne), order='F');

    # interior faces are zero
    f = zeros((nfe,ne),dtype=int,order='F');

    # face-to-element connectivity
    #f2t = mkf2t(t,elemtype,dim);
    f2t,t2t = mkf2e(t,elemtype,dim)[0:2];

    # find elements on the domain boundary
    ind = nonzero(f2t[2,:]==0)[0];

    for i in range(0,ind.size): # for each element on the domain boundary
        e = f2t[0,ind[i]]-1; # element e
        l = f2t[1,ind[i]]-1; # local face index
        for k in range(0,len(bndexpr)): # for each boundary expression
            a = bndexpr[k](reshape(pf[:,l,e],(dim,1))); # evaluate the boundary expression
            if a: # check if element e belong to this boundary
                f[l,e] = k+1; # then set f(l,e) to k
                break

    nprd = len(prdexpr);
    if nprd>0:
        f = f.flatten('F');
        pf = reshape(pf, (dim, nfe*ne), order = 'F');
        tprd = t.flatten('F');
        # periodic boundary faces are negative
        for i in range(0,nprd):
            i1 = nonzero(f == prdexpr[i][0])[0];
            f[i1] = -f[i1];            
            p1 = prdexpr[i][1](pf[:,i1]);
            e1 = int64(ceil(float64(i1+1)/float(nfe))) - 1;      # 1st elements
            l1 = i1 - e1*nfe;   # 1st local faces

            i2 = nonzero(f == prdexpr[i][2])[0];
            f[i2] = -f[i2];
            p2 = prdexpr[i][3](pf[:,i2]);
            e2 = int64(ceil(float64(i2+1)/float(nfe))) - 1;      # 2nd elements
            l2 = i2 - e2*nfe;   # 2nd local faces

            # update t2t to connect periodic elements
            for j in range(0,len(e1)):
                t2t[l2[j],e2[j]] = e1[j];
                t2t[l1[j],e1[j]] = e2[j];

            p1 = reshape(p1, (int(p1.size/len(i1)), len(i1)), 'F');
            p2 = reshape(p2, (int(p2.size/len(i2)), len(i2)), 'F');
            ind = xiny(p1,p2,0);

            v1 = t2fl[:,i1];
            v1 = sort(unique(v1.flatten()));
            p1 = prdexpr[i][1](p[:,v1]);

            v2 = t2fl[:,i2[ind]];
            v2 = sort(unique(v2.flatten()));
            p2 = prdexpr[i][3](p[:,v2]);

            p1 = reshape(p1, (int(p1.size/len(v1)), len(v1)), 'F');
            p2 = reshape(p2, (int(p2.size/len(v2)), len(v2)), 'F');
            ind = xiny(p1,p2,0);
            v2 = v2[ind];

            for j in range(0, v1.size):    
                i1 = nonzero(tprd == v2[j])[0];                           
                tprd[i1] = v1[j];
        tprd = reshape(tprd, (nve,ne), 'F');
    else:
        tprd = t;

    f = reshape(f, (nfe, ne), order = 'F');

    return f, tprd, t2t
