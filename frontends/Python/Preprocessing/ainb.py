from numpy import *

def ainb(x,y):

    m = len(x);
    ind = zeros(m);

    for j in range(0,m):
        d2 = abs(y-x[j])
        id = argmin(d2);
        if d2[id]<1e-8:
            ind[j]=id;

    return ind
