from numpy import *
from scipy.special import *

def koornwinder2d(x,n):

    x = 2.0*array(x)-1
    m = x.shape[0]
    k = int((n+1)*(n+2)/2)

    f = zeros((m,k))
    fx = zeros((m,k))
    fy = zeros((m,k))

    pq = pascalindex2d(k)

    xc = x.copy()
    xc[:,1] = minimum(0.99999999, xc[:,1])
    e = xc.copy()
    e[:,0] = 2*(1+xc[:,0])/(1-xc[:,1])-1
    i1 = where(x[:,1] == 1)
    e[i1,0] = -1
    e[i1,1] =  1

    for ii in range(0,k):
        pp = jacobi(pq[ii,0],0,0)
        qp = jacobi(pq[ii,1],2*pq[ii,0]+1,0)
        for i in range(0,int(pq[ii,0])):
            qp = convolve([-0.5,0.5],qp)

        pval = polyval(pp,e[:,0])
        qval = polyval(qp,e[:,1])

        fc = sqrt((2*pq[ii,0]+1)*2*(pq[ii,0]+pq[ii,1]+1))
        f[:,ii] = fc*(pval*qval)

    e[:,0] = 2*(1+xc[:,0])/(1-xc[:,1])-1
    e[:,1] = xc[:,1]
    de1 = xc.copy()
    de1[:,0] = 2./(1-xc[:,1])
    de1[:,1] = 2*(1+xc[:,0])/((1-xc[:,1])*(1-xc[:,1]))

    for ii in range(0,k):
        pp = jacobi(pq[ii,0],0,0)
        qp = jacobi(pq[ii,1],2*pq[ii,0]+1,0)
        for i in range(0,int(pq[ii,0])):
            qp = convolve([-0.5,0.5],qp)

        dpp = polyder(pp)
        dqp = polyder(qp)

        pval = polyval(pp,e[:,0])
        qval = polyval(qp,e[:,1])
        dpval = polyval(dpp,e[:,0])
        dqval = polyval(dqp,e[:,1])

        fc = sqrt((2*pq[ii,0]+1)*2*(pq[ii,0]+pq[ii,1]+1))
        fx[:,ii] = fc*dpval*qval*de1[:,0]
        fy[:,ii] = fc*(dpval*qval*de1[:,1] + pval*dqval)

    fx = 2*fx
    fy = 2*fy
    fz = 0

    return f, fx, fy, fz


def pascalindex2d(p):
    pq = zeros((p,2))
    l=0
    for i in range(0,p+1):
        for j in range(0,i+1):
            pq[l,0]=i-j
            pq[l,1]=j
            l = l+1
            if (l>=p):
                return pq
