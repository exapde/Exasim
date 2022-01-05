from numpy import array
from sympy import exp, sin, pi

def flux(u, q, w, v, x, t, mu, eta):
    f = mu[0]*q + u[0]*mu[1:];
    return f;

def source(u, q, w, v, x, t, mu, eta):
    x1 = x[0];
    x2 = x[1];
    v1 = mu[1];
    v2 = mu[2];
    s = array([(v1*x2*exp(v1)*exp(v2)*(exp(v2*(x2 - 1)) - 1)*(exp(v1*(x1 - 1)) + v1*x1*exp(v1*(x1 - 1)) - 1))/((exp(v1) - 1)*(exp(v2) - 1)) - (v2*x1*exp(v1 + v2*x2)*(exp(v1*(x1 - 1)) - 1)*(v2*x2 + 2))/((exp(v1) - 1)*(exp(v2) - 1)) - (v1*x2*exp(v2 + v1*x1)*(exp(v2*(x2 - 1)) - 1)*(v1*x1 + 2))/((exp(v1) - 1)*(exp(v2) - 1)) + (v2*x1*exp(v1)*exp(v2)*(exp(v1*(x1 - 1)) - 1)*(exp(v2*(x2 - 1)) + v2*x2*exp(v2*(x2 - 1)) - 1))/((exp(v1) - 1)*(exp(v2) - 1))]);
    return s;

def fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    f = flux(u, q, w, v, x, t, mu, eta);
    fb = array([f[0]*n[0] + f[1]*n[1] + tau[0]*(u[0]-uhat[0])]);
    return fb;

def ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    ub = array([0.0]);
    return ub;

def initu(x, mu, eta):
    u0 = array([0.0]);
    return u0;
