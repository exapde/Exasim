from numpy import array
from sympy import sin, pi, exp

def mass(u, q, w, v, x, t, mu, eta):
    m = array([1.0]);
    return m;

def flux(u, q, w, v, x, t, mu, eta):
    f = array([mu[0]*x[1]*u[0], -mu[1]*x[0]*u[0]]);
    return f;

def source(u, q, w, v, x, t, mu, eta):
    s = array([0.0]);
    return s;

def fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    f = flux(u, q, w, v, x, t, mu, eta);
    fb = array([f[0]*n[0] + f[1]*n[1] + tau[0]*(u[0]-0.0)]);
    return fb;

def ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    ub = array([0.0]);
    return ub;

def initu(x, mu, eta):
    u0 = array([exp(-(x[0]-0.25)*(x[0]-0.25)/0.01-(x[1]-0.0)*(x[1]-0.0)/0.01)]);
    return u0;
