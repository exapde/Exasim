from numpy import *

def xiny(x,y,opt):
# Determine if each row of x is a member of y
# If row j of x is a member of y and x(j,:) = y(k,:) then in[j] = k
# Else in[j] = 0

    if opt == 1:
        m = x.shape[0];
        dim = x.shape[1];
    else:
        dim = x.shape[0];
        m = x.shape[1];

    ind = zeros(m).astype(int);
    if opt==1:
        if dim==1:
            for j in range(0,m):
                d2 = (y[:,0]-x[j,0])**2
                id = d2.argmin();
                if d2[id]<1e-12:
                    ind[j]=id;
        elif dim==2:
            for j in range(0,m):
                d2 = (y[:,0]-x[j,0])**2 + (y[:,1]-x[j,1])**2
                id = d2.argmin()
                if d2[id]<1e-12:
                    ind[j] = id
        elif dim==3:
            for j in range(0,m):
                d2 = (y[:,0]-x[j,0])**2 + (y[:,1]-x[j,1])**2 + (y[:,2]-x[j,2])**2
                id = d2.argmin()
                if d2[id]<1e-12:
                    ind[j] = id
        else:
            n = y.shape[0]
            for j in range(0,m):
                d2 = sum((y - tile(x[j,:],(n,1)))**2, axis=1)
                id = d2.argmin()
                if d2[id]<1e-12:
                    ind[j] = id
    else:
        if dim==1:
            for j in range(0,m):
                d2 = (y[0,:]-x[0,j])**2
                id = d2.argmin();
                if d2[id]<1e-12:
                    ind[j]=id;
        elif dim==2:
            for j in range(0,m):
                d2 = (y[0,:]-x[0,j])**2 + (y[1,:]-x[1,j])**2
                id = d2.argmin()
                if d2[id]<1e-12:
                    ind[j] = id
        elif dim==3:
            for j in range(0,m):
                d2 = (y[0,:]-x[0,j])**2 + (y[1,:]-x[1,j])**2 + (y[2,:]-x[2,j])**2
                id = d2.argmin()
                if d2[id]<1e-12:
                    ind[j] = id
        else:
            n = y.shape[1]
            for j in range(0,m):
                d2 = sum((y - tile(x[:,j],(1,n)))**2, axis=0)
                id = d2.argmin()
                if d2[id]<1e-12:
                    ind[j] = id

    return ind
