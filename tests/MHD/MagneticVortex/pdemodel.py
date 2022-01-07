from numpy import array, reshape
from sympy import sin, pi, exp, sqrt

def mass(u, q, w, v, x, t, mu, eta):
    m = array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);
    return m;

def flux(u, q, w, v, x, t, mu, eta):
    gam = mu[0];
    gam1 = gam - 1.0;
    r = u[0];
    ru = u[1];
    rv = u[2];
    rE = u[3];
    bx = u[4];
    by = u[5];
    phi = u[6];
    r1 = 1/r;
    uv = ru*r1;
    vv = rv*r1;
    E = rE*r1;
    q = 0.5*(uv*uv+vv*vv);
    b = (bx*bx+by*by)*0.5;
    p = gam1*(rE-r*q-b-0.5*phi**2);
    h = E+p*r1+b*r1;
    uvb = uv*bx+vv*by;
    f = array([ru, ru*uv+p+b-bx*bx, rv*uv-by*bx, ru*h-uvb*bx+phi*bx, phi, -(vv*bx-by*uv), bx, rv, ru*vv-bx*by, rv*vv+p+b-by*by, rv*h-uvb*by+phi*by, -(uv*by-bx*vv), phi, by]);
    f = reshape(f,(7,2),'F');
    return f;

def source(u, q, w, v, x, t, mu, eta):
    s = array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    return s;

def fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    f = flux(u, q, w, v, x, t, mu, eta);
    fb = f[:,0]*n[0] + f[:,1]*n[1] + tau[0]*(u-uhat);
    return fb;

def ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    ub = array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    return ub;

def initu(x, mu, eta):
    x1 = x[0];
    x2 = x[1];
    gam = mu[0];
    r = sqrt(x1*x1+x2*x2);
    k = 1/(2*pi);    

    u01 = 1;
    u02 = u01-k*x2*exp(0.5*(1-r**2));
    u03 = u01+k*x1*exp(0.5*(1-r**2));
    u05 = -k*x2*exp(0.5*(1-r**2));
    u06 = k*x1*exp(0.5*(1-r**2));
    p = 1 + 1/(4*0.5)*(k**2*(1-2*0.5*r**2)-k**2*u01)*exp(2*0.5*(1-r**2));
    b = u06*u06 + u05*u05;
    v = u02*u02 + u03*u03;
    u04 = p/(gam-1) + u01/2*v + 0.5*b;
    u07 = 0;

    u0 = array([u01, u02, u03, u04, u05, u06, u07]);
    return u0;
