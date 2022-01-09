from numpy import *
from scipy.special import *

def koornwinder3d(x,n):

    x = 2.0*array(x)-1
    m = x.shape[0]
    k = int((n+1)*(n+2)*(n+3)/6)

    f = zeros((m,k))
    fx = zeros((m,k))
    fy = zeros((m,k))
    fz = zeros((m,k))

    pq = pascalindex3d(k)

    e = x.copy()
    for ii in range(0,m):
        if x[ii, 1] + x[ii, 2]:
            e[ii,0] = -2*(1+x[ii,0])/(x[ii,1]+x[ii,2])-1
        else:
            e[ii,0] = -1
        if x[ii,2] != 1:
            e[ii,1] = 2*(1+x[ii,1])/(1-x[ii,2])-1
        else:
            e[ii,1] = -1

    for ii in range(0,k):
        pp = jacobi(pq[ii,0],0,0)
        qp = jacobi(pq[ii,1],2*pq[ii,0]+1,0)
        rp = jacobi(pq[ii,2],2*pq[ii,0]+2*pq[ii,1]+2,0)

        for i in range(0,int(pq[ii,0])):
            qp = convolve([-0.5,0.5],qp)

        for i in range(0,int(pq[ii,0]+pq[ii,1])):
            rp = convolve([-0.5,0.5],rp)

        pval = polyval(pp,e[:,0])
        qval = polyval(qp,e[:,1])
        rval = polyval(rp,e[:,2])

        fc = sqrt((2*pq[ii,0]+1)*2*(pq[ii,0]+pq[ii,1]+1)*2*(pq[ii,0]+pq[ii,1]+pq[ii,2]+2))
        f[:,ii] = fc*(pval*qval*rval)

    xc = x.copy()
    for ii in range(0,m):
        if x[ii,1] + x[ii,2]==0:
            xc[ii,2] = -1e-8-xc[ii,1]
        if x[ii,2] == 1:
            xc[ii,2] = 0.99999999

    e[:,0] = -2.0*(1.0+xc[:,0])/(xc[:,1]+xc[:,2])-1.0
    e[:,1] = 2.0*(1.0+xc[:,1])/(1.0-xc[:,2])-1.0
    e[:,2] = xc[:,2]
    de1 = xc.copy()
    de1[:,0] = -2./(xc[:,1]+xc[:,2])
    de1[:,1] = 2*(1+xc[:,0])/((xc[:,1]+xc[:,2])*(xc[:,1]+xc[:,2]))
    de1[:,2] = 2*(1+xc[:,0])/((xc[:,1]+xc[:,2])*(xc[:,1]+xc[:,2]))
    de2 = xc.copy()
    de2[:,0] = 0*xc[:,0]
    de2[:,1] = 2/(1-xc[:,2])
    de2[:,2] = 2*(1+xc[:,1])/((1-xc[:,2])*(1-xc[:,2]))


    for ii in range(0,k):
        pp = jacobi(pq[ii,0],0,0)
        qp = jacobi(pq[ii,1],2*pq[ii,0]+1,0)
        rp = jacobi(pq[ii,2],2*pq[ii,0]+2*pq[ii,1]+2,0)

        for i in range(0,int(pq[ii,0])):
            qp = convolve([-0.5,0.5],qp)

        for i in range(0,int(pq[ii,0]+pq[ii,1])):
            rp = convolve([-0.5,0.5],rp)

        dpp = polyder(pp)
        dqp = polyder(qp)
        drp = polyder(rp)

        pval = polyval(pp,e[:,0])
        qval = polyval(qp,e[:,1])
        rval = polyval(rp,e[:,2])
        dpval = polyval(dpp,e[:,0])
        dqval = polyval(dqp,e[:,1])
        drval = polyval(drp,e[:,2])

        fc = sqrt((2*pq[ii,0]+1)*2*(pq[ii,0]+pq[ii,1]+1)*2*(pq[ii,0]+pq[ii,1]+pq[ii,2]+2))
        fx[:,ii] = fc*(dpval*qval*rval*de1[:,0])
        fy[:,ii] = fc*(dpval*qval*rval*de1[:,1] + pval*dqval*rval*de2[:,1])
        fz[:,ii] = fc*(dpval*qval*rval*de1[:,2] + pval*dqval*rval*de2[:,2] + pval*qval*drval)

    fx = 2.0*fx
    fy = 2.0*fy
    fz = 2.0*fz

    return f, fx, fy, fz


def pascalindex3d(p):
    pq = zeros((p,3))
    l=0
    for i in range(0,p+1):
        for j in range(0,i+1):
            for k in range(0,j+1):
                pq[l,0]=i-j
                pq[l,1]=j-k
                pq[l,2]=k
                l = l+1
                if (l>=p):
                    return pq
