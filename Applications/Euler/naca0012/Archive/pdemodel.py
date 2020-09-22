from numpy import array, reshape
from sympy import sin, pi, exp, sqrt

def mass(u, q, w, v, x, t, mu, eta):
    m = array([1.0, 1.0, 1.0, 1.0]);
    return m;

def flux(u, q, w, v, x, t, mu, eta):
    gam = mu[0];
    gam1 = gam - 1.0;
    r = u[0];
    ru = u[1];
    rv = u[2];
    rE = u[3];
    r1 = 1/r;
    uv = ru*r1;
    vv = rv*r1;
    E = rE*r1;
    p = gam1*(rE-r*0.5*(uv*uv+vv*vv));
    h = E+p*r1;
    f = array([ru, ru*uv+p, rv*uv, ru*h, rv, ru*vv, rv*vv+p, rv*h]);
    f = reshape(f,(4,2),'F');
    return f;

def source(u, q, w, v, x, t, mu, eta):
    s = array([0.0, 0.0, 0.0, 0.0]);
    return s;

def fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    f = flux(u, q, w, v, x, t, mu, eta);
    fb = f[:,0]*n[0] + f[:,1]*n[1] + tau[0]*(u-uhat);
    return fb;

def ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    ub = array([0.0, 0.0, 0.0, 0.0]);
    return ub;

def initu(x, mu, eta):
    t = 0.0;
    phi = 5.0;
    x1 = x[0];
    x2 = x[1];
    gam = mu[0];
    M_ref = mu[1];

    r = sqrt((x1-t)**2 + x2**2);
    u01 = (1 - ((gam-1)/(16*pi**2)) * phi**2 * exp(2*(1-r**2)) )**(1/(gam-1));
    u02 = u01 * (1 - M_ref**(-1.0) * phi * exp(1-r**2) * x2/(2*pi));
    u03 = u01 * M_ref**(-1.0) * phi * exp(1-r**2) * x1/(2*pi);
    p = u01**gam/(gam*M_ref**2);
    u04 = p/(gam-1) + 0.5 * (u02*u02 + u03*u03) / u01;

    u0 = array([u01, u02, u03, u04]);
    return u0;
