from numpy import *

def cubemesh(m,n,o,elemtype):
# Generate mesh for unit cube

    p = [[x/m, y/n, z/o] for z in arange(0,o+1) for y in arange(0,n+1) for x in arange(0,m+1)];
    p = array(p,order='F').T;

    m = m+1;
    n = n+1;
    o = o+1;

    if elemtype==0:
        t0 = [[5, 6, 7, 3], [5, 6, 2, 3], [5, 1, 2, 3], [6, 8, 7, 4], [6, 7, 4, 3], [6, 2, 4, 3]];
        t0 = array(t0);

        map = array([1, 2]);
        map = concatenate([map,  map+m]);
        map = concatenate([map, map+m*n]);
        t=map[t0-1];
        t=kron(t,ones((m-1,1)))+kron(ones(t.shape),reshape(arange(0,m-1),(m-1,1)));
        t=kron(t,ones((n-1,1)))+kron(ones(t.shape),reshape(arange(0,n-1)*m,(n-1,1)));
        t=kron(t,ones((o-1,1)))+kron(ones(t.shape),reshape(arange(0,o-1)*(m*n),(o-1,1)));
    else:
        t = array([1, 2, m+2, m+1, m*n+1, m*n+2, m*n+m+2, m*n+m+1]);
        t = kron(t,ones((o-1,1)))+kron(ones(t.shape),reshape(arange(0,o-1)*(m*n),(o-1,1)));
        t = kron(t,ones((n-1,1)))+kron(ones(t.shape),reshape(arange(0,n-1)*m,(n-1,1)));
        t = kron(t,ones((m-1,1)))+kron(ones(t.shape),reshape(arange(0,m-1),(m-1,1)));

    t = t.transpose().astype(int)-1;

    return p, t
