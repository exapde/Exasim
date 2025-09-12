from numpy import *
from getelemface import *
#import sys

def mkf2t(t,elemtype,dim):

    ne = t.shape[1]
    face = getelemface(dim,elemtype)
    nvf = face.shape[0]
    nfe = face.shape[1]

    N = nfe*ne;
    sfaces = reshape(t[face.flatten(order='F'),:], (nvf, N), order='F');
    sfaces = sort(sfaces, axis=0);
    f = 0*sfaces;
    f[:,0] = sfaces[:,0];
    f2t = zeros((4,N));
    f2t[0,0] = 1;
    f2t[1,0] = 1;

    k = 0;
    for i in range(1,N):
        l = (i+1) % nfe; # local face l
        if (l == 0):
            l = nfe;
        e = (i+1-l)/nfe + 1; # element e

        s = reshape(sfaces[:,i],(nvf,1));
        diff = sum(abs(s - f[:,0:(k+1)]),axis=0);
        dmin = min(diff);
        if dmin==0: # match
            j = nonzero(diff.flatten() == dmin)
            f2t[2,j] = e;
            f2t[3,j] = l;
            # print(diff.flatten())
            # print([k,i,j,e,l])
            # print(f[:,0:(k+1)])
            # sys.exit("here")
        else: # not mach
            k = k + 1;
            # add local face to the list
            f[:,k] = sfaces[:,i];
            f2t[0,k] = e;
            f2t[1,k] = l;

    f2t = f2t[:,0:(k+1)];
    f2t = f2t.astype(int);

    return f2t
