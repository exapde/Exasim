import numpy as np

def TriSquareMesh(m, n):

    p = [[x/m, y/n] for y in np.arange(0,n+1) for x in np.arange(0,m+1)];
    p = np.array(p,order='F');

    t = [ [ [i + (j - 1) * (m + 1), i + (j - 1) * (m + 1) + 1, i + j * (m + 1)], [i + (j - 1) * (m + 1) + 1, i + j * (m + 1) + 1, i + j * (m + 1)]] for i in range (1,m+1) for j in range (1,n+1) ]
    t = np.array(t);
    t = np.reshape(t,(2*m*n,3));

    return p, t

def QuadSquareMesh(m,n):

    p = [[x/m, y/n] for y in np.arange(0,n+1) for x in np.arange(0,m+1)];
    p = np.array(p,order='F');

    t = [ [i + (j - 1) * (m + 1), i + (j - 1) * (m + 1) + 1,  i + j * (m + 1) + 1,   i + j * (m + 1)] for j in range (1,m+1) for i in range (1,n+1)];
    t = np.array(t);

    return p, t

def SquareMesh(m, n, elemtype):

    if (elemtype==0):
        p, t = TriSquareMesh(m,n)[0:2]
    else:
        p, t = QuadSquareMesh(m,n)[0:2]

    p = p.transpose();
    t = t.transpose() - 1;

    return p, t
