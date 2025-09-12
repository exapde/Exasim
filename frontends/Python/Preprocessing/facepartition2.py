from numpy import *
from mkf2e import mkf2e
from faceconnectivity2 import faceconnectivity2

def facepartition2(dmd,t,f,bcm,dim,elemtype,porder,nproc):

    for i in range(0,nproc):
        fi = f[:,dmd[i]['elempart'].flatten()];
        f2t = mkf2e(t[:,dmd[i]['elempart'].flatten()],elemtype,dim)[0];

        # only on [intelem bndelem extelem]
        nelem = dmd[i]['elempartpts'];
        if (size(nelem)>=3):
            ne3 = sum(nelem[0:3]);
            ind = nonzero(f2t[0,:]<=ne3)[0];
            f2t = f2t[:,ind];

        # reorder so that boundary faces are last
        ina = nonzero(f2t[2,:]>0)[0]; # interior faces
        inb = nonzero(f2t[2,:]==0)[0]; # boundary faces
        #inc = sub2ind(shape(fi), f2t[2,inb], f2t[1,inb]);
        inc = (f2t[0,inb]-1)*(fi.shape[0]) + f2t[1,inb] - 1;
        fb = fi[unravel_index(inc, fi.shape, 'F')];

        fa = sort(unique(fb))-1; # boundary indices
        bcn = sort(unique(bcm[fa])); # a list of boundary conditions
        nbc = bcn.size;
        dmd[i]['facepartpts'] = array([ina.size]).astype(int);
        dmd[i]['facepartbnd'] = array([0]).astype(int);
        ind = zeros(fb.size).astype(int);
        m = 0;
        for j in range(0,nbc): # for each boundary condition bcn(j)
            bj = nonzero(bcm==bcn[j])[0] + 1; # find all boundaries that have condition bcn(j)
            n = 0;
            for k in range(0,bj.size): # for each boundary that has condition bcn(j)
                ii = nonzero(fb == bj[k])[0]; # indices of the boundary bj(k)
                l = ii.size;
                if l>0:
                    n = n + l;
                    ind[m:(m+l)] = ii;
                    m = m + l;
            dmd[i]['facepartpts'] = append(dmd[i]['facepartpts'], n); #concatenate([dmd[i]['facepartpts'].flatten(), [n]]);
            dmd[i]['facepartbnd'] = append(dmd[i]['facepartbnd'], bcn[j]);  #concatenate([dmd[i]['facepartbnd'].flatten(), [bcn[j]]]);

        # [interior faces, boundary faces]
        f2t = f2t[:,concatenate([ina, inb[ind]])];
        dmd[i]['facepartpts'] = array(dmd[i]['facepartpts']);
        dmd[i]['facepartbnd'] = array(dmd[i]['facepartbnd']);
        dmd[i]['facecon'] = faceconnectivity2(t[:,dmd[i]['elempart'].flatten('F')],f2t,dim,elemtype,porder)[0];
        
    return dmd
