from numpy import array
from sympy import exp, sin, pi

def flux(u, q, w, v, x, t, mu, eta):
    f = mu[0]*q + u[0]*mu[1:3] + u[0]*u[0]*mu[3:];
    return f;

def source(u, q, w, v, x, t, mu, eta):
    s = array([0.0]);
    return s;

def fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    f = flux(u, q, w, v, x, t, mu, eta);
    tm = f[0]*n[0] + f[1]*n[1] + tau[0]*(u[0]-uhat[0]);
    fb = array([tm, tm]);
    return fb;

def ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    ub = array([1-2*x[0], u]);
    return ub;

def initu(x, mu, eta):
    u0 = array([1-2*x[0]]);
    return u0;
