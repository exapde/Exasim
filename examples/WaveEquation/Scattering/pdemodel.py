from numpy import array
from sympy import sin, cos, sqrt, pi

def mass(u, q, w, v, x, t, mu, eta):
    m = array([1.0]);
    return m;

def flux(u, q, w, v, x, t, mu, eta):
    f = mu[0]*q;
    return f;

def source(u, q, w, v, x, t, mu, eta):
    s = array([0.0]);
    return s;

def fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    #ui = sin(mu[1]*x[0]+mu[2]*x[1]-mu[3]*sqrt(mu[1]*mu[1]+mu[2]*mu[2])*t);
    uix = mu[1]*cos(mu[1]*x[0]+mu[2]*x[1]-mu[3]*sqrt(mu[1]*mu[1]+mu[2]*mu[2])*t);
    uiy = mu[2]*cos(mu[1]*x[0]+mu[2]*x[1]-mu[3]*sqrt(mu[1]*mu[1]+mu[2]*mu[2])*t);
    fb = array([u, uix*n[0]+uiy*n[1]]); # Neumann conditions
    return fb;

def ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    ub = array([u, u]);
    return ub;

def initu(x, mu, eta):
    u0 = array([0.0]);
    return u0;

def initq(x, mu, eta):
    q0 = array([0.0, 0.0]);
    return q0;

def initw(x, mu, eta):
    w0 = array([0.0]);
    return w0;
