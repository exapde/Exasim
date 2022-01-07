from numpy import array, reshape, copy, hstack, dot
from sympy import sin, pi, exp, sqrt, tanh

def mass(u, q, w, v, x, t, mu, eta):
    m = array([1.0, 1.0, 1.0, 1.0]);
    return m;

def flux(u, q, w, v, x, t, mu, eta):
    gam = mu[0];
    gam1 = gam - 1.0;
    Re = mu[1];
    Pr = mu[2];
    Minf = mu[3];
    Re1 = 1/Re;
    M2 = Minf**2;
    c23 = 2.0/3.0;
    fc = 1/(gam1*M2*Re*Pr);
    r = u[0];
    ru = u[1];
    rv = u[2];
    rE = u[3];
    rx = q[0];
    rux = q[1];
    rvx = q[2];
    rEx = q[3];
    ry = q[4];
    ruy = q[5];
    rvy = q[6];
    rEy = q[7];
    r1 = 1/r;
    uv = ru*r1;
    vv = rv*r1;
    E = rE*r1;
    ke = 0.5*(uv*uv+vv*vv);
    p = gam1*(rE-r*ke);
    h = E+p*r1;
    fi = array([ru, ru*uv+p, rv*uv, ru*h, rv, ru*vv, rv*vv+p, rv*h]);
    ux = (rux - rx*uv)*r1;
    vx = (rvx - rx*vv)*r1;
    kex = uv*ux + vv*vx;
    px = gam1*(rEx - rx*ke - r*kex);
    Tx = gam*M2*(px*r - p*rx)*r1**2;
    uy = (ruy - ry*uv)*r1;
    vy = (rvy - ry*vv)*r1;
    key = uv*uy + vv*vy;
    py = gam1*(rEy - ry*ke - r*key);
    Ty = gam*M2*(py*r - p*ry)*r1**2;
    txx = Re1*c23*(2*ux - vy);
    txy = Re1*(uy + vx);
    tyy = Re1*c23*(2*vy - ux);
    fv = array([0, txx, txy, uv*txx + vv*txy + fc*Tx, 0, txy, tyy, uv*txy + vv*tyy + fc*Ty]);
    f = fi+fv;
    f = reshape(f,(4,2),'F');
    return f;

def source(u, q, w, v, x, t, mu, eta):
    s = array([0.0, 0.0, 0.0, 0.0]);
    return s;

def fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    f = flux(u, q, w, v, x, t, mu, eta);
    fi = f[:,0]*n[0] + f[:,1]*n[1] + tau[0]*(u-uhat); # numerical flux at freestream boundary
    fw = fi.flatten('F');
    fw[0] = 0.0; # zero velocity
    fw[3] = 0.0; # adiabatic wall -> zero heat flux
    fb = hstack((fw, fi)); # wall and freestream boundary conditions
    fb = reshape(fb,(4,2),'F');
    return fb;

def ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
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
    uinf = mu[4:8]; # freestream flow
    uinf = uinf[:];
    u = u[:];          # state variables
    ui = 0.5*((u+uinf) + dot(An,(u-uinf)));  # Riemann solution
    uw = u.flatten('F');
    uw[1] = 0; # zero velocity
    uw[2] = 0; # zero velocity

    ub = hstack((uw, ui)); # wall and freestream boundary conditions
    ub = reshape(ub,(4,2),'F');
    return ub;

def initu(x, mu, eta):
    u01 = mu[4];
    u02 = mu[5];
    u03 = mu[6];
    u04 = mu[7];

    u0 = array([u01, u02, u03, u04]);
    return u0;
