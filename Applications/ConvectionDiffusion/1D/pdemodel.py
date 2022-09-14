from numpy import array
from sympy import sin, cos

def flux(u, q, w, v, x, t, mu, eta):
    f = mu[0]*q + u[0]*mu[1]
    return f

def source(u, q, w, v, x, t, mu, eta):
    x1 = x[0]
    nu = mu[0]
    c = mu[1]
    c2 = c**2

    s = array([c2*nu*(sin(c*x1) - cos(c*x1))])
    return s

def fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    f = flux(u, q, w, v, x, t, mu, eta)
    fb = array(f[0]*n[0] + tau*(u[0]-uhat[0]))
    return fb

def ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    ub = array([0.0])
    return ub

def initu(x, mu, eta):
    u0 = array([0.0])
    return u0