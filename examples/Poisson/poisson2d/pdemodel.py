
from numpy import array
from sympy import sin, pi

def flux(u, q, w, v, x, t, mu, eta):
    f = mu*q;
    return f;

def source(u, q, w, v, x, t, mu, eta):
    x1 = x[0];
    x2 = x[1];
    s = array([(2*pi*pi)*sin(pi*x1)*sin(pi*x2)]);
    return s;

def visscalars(u, q, w, v, x, t, mu, eta):
    s = 0*q;
    s[0] = u[0];
    s[1] = q[0] + q[1];
    return s;

def visvectors(u, q, w, v, x, t, mu, eta):    
    return q;

def qoivolume(u, q, w, v, x, t, mu, eta):    
    s = 0*q;
    x1 = x[0];
    x2 = x[1];
    uexact = sin(pi*x1)*sin(pi*x2);
    s[0] = (u[0] - uexact)*(u[0] - uexact);
    s[1] = u[0];
    return s;    

def fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    f = flux(u, q, w, v, x, t, mu, eta);
    fb = array([f[0]*n[0] + f[1]*n[1] + tau[0]*(u[0]-0.0)]);
    return fb;

def qoiboundary(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    f = flux(u, q, w, v, x, t, mu, eta);
    fb = array([f[0]*n[0] + f[1]*n[1] + tau[0]*(u[0]-uhat[0])]);
    return fb;

def fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    fb = array([tau[0]*(0.0-uhat[0])]);
    return fb;

def ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    ub = array([0.0]);
    return ub;

def initu(x, mu, eta):
    u0 = array([0.0]);
    return u0;
