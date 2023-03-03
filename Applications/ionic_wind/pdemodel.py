from numpy import array
from sympy import exp, sin, pi

def mass(u, q, w, v, x, t, mu, eta):
    m = array([1.0]);
    return m;

def flux(u, q, w, v, x, t, mu, eta):
    r = x[0]
    z = x[1]
    f = mu[0]*q + u[0]*mu[1:3];      # mu1*q + c dot u
    return f;

def source(u, q, w, v, x, t, mu, eta):
    s = array([0.0])
    return s;

def fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    f = flux(u, q, w, v, x, t, mu, eta);
    f_dot_n = f[0]*n[0] + f[1]*n[1] + tau[0]*(u[0]-uhat[0])
    fb = array([f_dot_n]);
    return fb;

def ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    ub = array([0.0]);      # Dirichlet BC on all sides
    return ub;

def initu(x, mu, eta):
    x1 = x[0]
    y1 = x[1]
    sigx = 0.125
    sigy = 0.125
    x0 = 0.5
    y0 = 0.5

    u0 = array([exp(-0.5*( (x1-x0)**2/sigx**2 + (y1-y0)**2/sigy**2) )]);
    return u0;
