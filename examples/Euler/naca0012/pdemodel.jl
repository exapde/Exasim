function mass(u, q, w, v, x, t, mu, eta)
    m = [1.0, 1.0, 1.0, 1.0];
    return m;
end
function flux(u, q, w, v, x, t, mu, eta)
    gam = mu[1];
    gam1 = gam - 1.0;
    r = u[1];
    ru = u[2];
    rv = u[3];
    rE = u[4];
    r1 = 1/r;
    uv = ru*r1;
    vv = rv*r1;
    E = rE*r1;
    p = gam1*(rE-r*0.5*(uv*uv+vv*vv));
    h = E+p*r1;
    f = [ru, ru*uv+p, rv*uv, ru*h, rv, ru*vv, rv*vv+p, rv*h];
    f = reshape(f,(4,2));
    return f;
end
function source(u, q, w, v, x, t, mu, eta)
    s = [0.0, 0.0, 0.0, 0.0];
    return s;
end
function ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    ub = zeros((4,2));
    return ub;
end
function fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    f = flux(u, q, w, v, x, t, mu, eta);

    gam = mu[1];
    gam1 = gam - 1.0;
    r = u[1];
    ru = u[2];
    rv = u[3];
    rE = u[4];
    nx = n[1];
    ny = n[2];

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

    K = [ 1  1  0  1 ; un-a  un  0  un+a ; ut  ut  1  ut ; h - un*a  (1/2)*(un^2 + ut^2)  ut  h+un*a];
    Kinv = (gam1/(2*a^2))*[ h + (a/gam1)*(un-a)  -(un+a/gam1)  -ut  1; -2*h + (4/gam1)*a^2  2*un  2*ut  -2 ; -2*(ut*a^2)/gam1  0  2*(a^2)/gam1  0 ; h - a*(un+a)/gam1  -un+a/gam1  -ut  1];
    T = [ 1 0 0 0; 0 nx ny 0; 0 -ny nx 0; 0 0 0 1];
    Tinv = [ 1 0 0 0; 0 nx -ny 0; 0 ny nx 0; 0 0 0 1];
    Lambda = [ tanh(1e3*(un-a)) 0 0 0 ; 0 tanh(1e3*(un)) 0 0 ; 0 0 tanh(1e3*(un)) 0 ; 0 0 0 tanh(1e3*(un+a)) ];
    # E = simplify(K * Lambda * Kinv);
    # An = simplify(Tinv * E * T);
    E = (K * Lambda * Kinv);
    An = (Tinv * E * T);

    # freestream boundary condition
    uinf = mu[3:6]; # freestream flow
    uinf = uinf[:];
    u = u[:];          # state variables
    ui = 0.5*((u.+uinf) + An*(u.-uinf));  # Riemann solution
    fi = f[:,1]*n[1] + f[:,2]*n[2] + tau[1]*(u-ui);

    uw = u[:];
    uw[2] = ru - run*nx;
    uw[3] = rv - run*ny;
    uw = 0.5*((u.+uw) + An*(u.-uw));  # Riemann solution
    fw = f[:,1]*n[1] + f[:,2]*n[2] + tau[1]*(u-uw); # numerical flux at freestream boundary
    fw[1] = 0;
    fw[end] = 0;
    fb = hcat(fw[:], fi[:]); # wall and freestream boundary conditions

    return fb;
end
function initu(x, mu, eta)

    u0 = mu[3:6]; # freestream flow
    return u0;
end
