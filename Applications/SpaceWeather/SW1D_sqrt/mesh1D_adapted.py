from numpy import *
from scipy.optimize import fsolve

def mesh1D_adapted(r1,r2):
    dlay = 1
    dwall = 1e-2
    nx = 16

    #Calculate mesh ratio
    c = 1-dlay/dwall
    def function(x):
        return scalingfun(x,nx,c)

    rat = fsolve(function,1.1)

    xv = zeros(nx+1)
    xv[1] = dwall
    for i in range(0,nx-1):
        xv[i+2] = xv[i+1] + dwall*rat**i

    #if abs(xv[-1]-dlay)>1e-8:
        #error("Something wrong with the input parameters")

    p = array([xv*(r2-r1)+r1])
    t = array([[k,k+1] for k in range(0,nx)])
    t = t.transpose()

    return p,t


def scalingfun(x,n,c):
    F = c

    for i in range(0,n-1):
        F += x**i 

    return F