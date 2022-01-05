from numpy import array
from sympy import sin, pi

def flux(u, q, w, v, x, t, mu, eta):
    f = mu[0]*q;
    return f;

def source(u, q, w, v, x, t, mu, eta):
    s = array([0.0]);
    return s;

def fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    f = flux(u, q, w, v, x, t, mu, eta);
    tm = array([f[0]*n[0] + f[1]*n[1] + f[2]*n[2]+ tau[0]*(u[0]-uhat[0])]);
    fb = array([tm[0], mu[2]*n[0]+mu[3]*n[1]+mu[4]*n[2]]);
    return fb;

def ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    ub = array([mu[1], u]);
    return ub;

def initu(x, mu, eta):
    u0 = array([1.0]);
    return u0;
