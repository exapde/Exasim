from numpy import *
from scipy.special import *

def koornwinder1d(x,n):

    x = 2.0*array(x)-1
    np = len(x)
    x = reshape(x,(np,1))
    f = zeros((np,n+1))
    fx = zeros((np,n+1))
    for ii in range(0,n+1):
        P = jacobi(ii,0,0)
        P = P*sqrt(2*ii+1)
        D = polyder(P)
        f[:,ii] =  polyval(P,x)[:,0]
        fx[:,ii] = polyval(D,x)[:,0]

    fx = 2*fx
    fy = 0
    fz = 0

    return f, fx, fy, fz
