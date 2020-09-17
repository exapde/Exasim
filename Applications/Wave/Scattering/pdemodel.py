from numpy import array
from sympy import sin, cos, pi, exp

def flux(u, q, w, v, x, t, mu, eta):
    z = x[0];
    r = x[1];
    f = mu*r*q;
    return f;

def source(u, q, w, v, x, t, mu, eta):
    z = x[0];
    r = x[1];
    s = array([sin(r)/exp(z)]);
    return s;

def fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    f = flux(u, q, w, v, x, t, mu, eta);
    tm = f[0]*n[0] + f[1]*n[1] + tau[0]*(u[0]-uhat[0]);
    fb = array([0.0, tm, tm, tm]);
    return fb;

def ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    z = x[0];
    r = x[1];
    uexact = exp(-z)*cos(r);
    ub = array([u, uexact, uexact, uexact]);     
    return ub;

def initu(x, mu, eta):
    u0 = array([0.0]);
    return u0;
