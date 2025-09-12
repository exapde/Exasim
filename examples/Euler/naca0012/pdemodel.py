from numpy import array, reshape, zeros, dot, hstack
from sympy import sin, pi, exp, sqrt, tanh

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

    gam = mu[0];
    gam1 = gam - 1.0;
    r = u[0];
    ru = u[1];
    rv = u[2];
    rE = u[3];
    nx = n[0];
    ny = n[1];

    r1 = 1/r;
    uv = ru*r1;
    vv = rv*r1;
    E = rE*r1;
    p = gam1*(rE-r*0.5*(uv*uv+vv*vv));
    h = E+p*r1;
    a = sqrt(gam*p*r1);

    run = ru*nx + rv*ny;
    rut = -ru*ny + rv*nx;
    un = run/r;
    ut = rut/r;

    K = array([[ 1,  1,  0,  1], [un-a,  un,  0,  un+a ], [ut,  ut,  1,  ut], [h-un*a,  (1/2)*(un**2 + ut**2),  ut,  h+un*a]]);
    Kinv = (gam1/(2*a**2))*array([[h+(a/gam1)*(un-a),  -(un+a/gam1),  -ut,  1], [-2*h+(4/gam1)*a**2,  2*un,  2*ut,  -2], [-2*(ut*a**2)/gam1,  0,  2*(a**2)/gam1,  0], [h-a*(un+a)/gam1,  -un+a/gam1,  -ut,  1]]);
    T = array([[ 1, 0, 0, 0], [0, nx, ny, 0], [0, -ny, nx, 0], [0, 0, 0, 1]]);
    Tinv = array([[ 1, 0, 0, 0], [0, nx, -ny, 0], [0, ny, nx, 0], [0, 0, 0, 1]]);
    Lambda = array([[tanh(1e3*(un-a)), 0, 0, 0], [0, tanh(1e3*(un)), 0, 0 ], [0, 0, tanh(1e3*(un)), 0 ], [0, 0, 0, tanh(1e3*(un+a)) ]]);
    E = (K * Lambda * Kinv);
    An = (Tinv * E * T);

    # freestream boundary condition
    uinf = mu[2:6]; # freestream flow
    uinf = uinf[:];
    u = u[:];          # state variables
    ui = 0.5*((u+uinf) + dot(An,(u-uinf)));  # Riemann solution
    fi = f[:,0]*n[0] + f[:,1]*n[1] + tau[0]*(u-ui);

    uw = u.flatten('F');
    uw[1] = ru - run*nx;
    uw[2] = rv - run*ny;
    uw = 0.5*((u+uw) + dot(An,(u-uw)));  # Riemann solution
    fw = f[:,0]*n[0] + f[:,1]*n[1] + tau[0]*(u-uw);
    fw[0] = 0;
    fw[3] = 0;

    fb = hstack((fw, fi)); # wall and freestream boundary conditions
    fb = reshape(fb,(4,2),'F');
    return fb;

def ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    ub = zeros((4,2));
    return ub;

def initu(x, mu, eta):
    u01 = mu[2];
    u02 = mu[3];
    u03 = mu[4];
    u04 = mu[5];

    u0 = array([u01, u02, u03, u04]);
    return u0;
