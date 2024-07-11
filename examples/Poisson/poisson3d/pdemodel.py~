from numpy import array
from sympy import sin, pi

def flux(u, q, w, v, x, t, mu, eta):
    f = mu[0]*q;
    return f;

def source(u, q, w, v, x, t, mu, eta):
    x1 = x[0];
    x2 = x[1];
    x3 = x[2];
    s = array([(3*pi*pi)*sin(pi*x1)*sin(pi*x2)*sin(pi*x3)]);
    return s;

def fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    f = flux(u, q, w, v, x, t, mu, eta);
    fb = array([f[0]*n[0] + f[1]*n[1] + f[2]*n[2]+ tau[0]*(u[0]-0.0)]);
    return fb;

def fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau):    
    fb = array([tau[0]*(0.0-uhat[0])]);
    return fb;

def ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    ub = array([mu[1]]);
    return ub;

def initu(x, mu, eta):
    u0 = array([0.0]);
    return u0;
