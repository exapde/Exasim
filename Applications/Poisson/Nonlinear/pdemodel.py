
from numpy import array
from sympy import sin, pi

def flux(u, q, w, v, x, t, mu, eta):
    f = (1+u*u)*mu[0]*q;
    return f;

def source(u, q, w, v, x, t, mu, eta):
    x1 = x[0];
    x2 = x[1];
    s = array([x1*sin(x2)]);
    return s;

def fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    f = flux(u, q, w, v, x, t, mu, eta);
    tm = f[0]*n[0] + f[1]*n[1] + tau[0]*(u[0]-uhat[0]);
    fb = array([tm, 0.0]);    
    return fb;

def ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    ub = array([mu[1], u]);
    return ub;

def initu(x, mu, eta):
    u0 = array([1.0]);
    return u0;
