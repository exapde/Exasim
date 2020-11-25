
function pde = pdemodel
pde.mass = @mass;
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.ubou = @ubou;
pde.initu = @initu;
end

function m = mass(u, q, w, v, x, t, mu, eta)
m = sym([1.0; 1.0; 1.0; 1.0; 1.0]);
end

function f = flux(u, q, w, v, x, t, mu, eta)
    gam = mu(1);
    gam1 = gam - 1.0;
    r = u(1);
    ru = u(2);
    rv = u(3);
    rw = u(4);
    rE = u(5);
    r1 = 1/r;
    uv = ru*r1;
    vv = rv*r1;
    wv = rw*r1;
    E = rE*r1;
    p = gam1*(rE-r*0.5*(uv*uv+vv*vv+wv*wv));
    h = E+p*r1;
    f = [ru, ru*uv+p, rv*uv, rw*uv, ru*h, ...
          rv, ru*vv, rv*vv+p, rw*vv, rv*h, ...
          rw, ru*wv, rv*wv, rw*wv+p, rw*h];
    f = reshape(f,[5,3]);  
end
function f = flux(u, q, w, v, x, t, mu, eta)
    gam = mu(1);    
    gam1 = gam - 1.0;
    Re = mu(2);
    Pr = mu(3);
    Minf = mu(4);
    Re1 = 1/Re;
    M2 = Minf^2;
    c23 = 2.0/3.0;
    fc = 1/(gam1*M2*Re*Pr);
    r = u(1);
    ru = u(2);
    rv = u(3);
    rw = u(4);
    rE = u(5);
    rx = q(1);
    rux = q(2);
    rvx = q(3);
    rwx = q(4);
    rEx = q(5);
    ry = q(6);
    ruy = q(7);
    rvy = q(8);
    rwy = q(9);
    rEy = q(10);
    rz = q(11);
    ruz = q(12);
    rvz = q(13);
    rwz = q(14);
    rEz = q(15);
    r1 = 1/r;
    uv = ru*r1;
    vv = rv*r1;
    wv = rw*r1;
    E = rE*r1;
    ke = 0.5*(uv*uv+vv*vv+wv*wv);
    p = gam1*(rE-r*ke);
    h = E+p*r1;
    fi = [ru, ru*uv+p, rv*uv, rw*uv, ru*h, ...
          rv, ru*vv, rv*vv+p, rw*vv, rv*h, ...
          rw, ru*wv, rv*wv, rw*wv+p, rw*h];
    ux = (rux - rx*uv)*r1;
    vx = (rvx - rx*vv)*r1;
    wx = (rwx - rx*wv)*r1;
    kex = uv*ux + vv*vx + wv*wx;
    px = gam1*(rEx - rx*ke - r*kex);
    Tx = gam*M2*(px*r - p*rx)*r1^2;
    uy = (ruy - ry*uv)*r1;
    vy = (rvy - ry*vv)*r1;
    wy = (rwy - ry*wv)*r1;
    key = uv*uy + vv*vy + wv*wy;
    py = gam1*(rEy - ry*ke - r*key);
    Ty = gam*M2*(py*r - p*ry)*r1^2;
    uz = (ruz - rz*uv)*r1;
    vz = (rvz - rz*vv)*r1;
    wz = (rwz - rz*wv)*r1;
    kez = uv*uz + vv*vz + wv*wz;
    pz = gam1*(rEz - rz*ke - r*kez);
    Tz = gam*M2*(pz*r - p*rz)*r1^2;
    txx = Re1*c23*(2*ux - vy - wz);
    txy = Re1*(uy + vx);
    txz = Re1*(uz + wx);
    tyy = Re1*c23*(2*vy - ux - wz);
    tyz = Re1*(vz + wy);
    tzz = Re1*c23*(2*wz - ux - vy);
    fv = [0, txx, txy, txz, uv*txx + vv*txy + wv*txz + fc*Tx, ...
          0, txy, tyy, tyz, uv*txy + vv*tyy + wv*tyz + fc*Ty, ...
          0, txz, tyz, tzz, uv*txz + vv*tyz + wv*tzz + fc*Tz];
    f = fi+fv;
    f = reshape(f,[5,3]);    
end


function s = source(u, q, w, v, x, t, mu, eta)
    s = [sym(0.0); sym(0.0); sym(0.0); sym(0.0); sym(0.0)];
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    f = flux(u, q, w, v, x, t, mu, eta);
    
    gam = mu(1);
    gam1 = gam - 1.0;
    
    nx = n(1);
    ny = n(2);
    nz = n(3);
    
    r = u(1);
    ru = u(2);
    rv = u(3);
    rw = u(4);
    rE = u(5);
    
    r1 = 1/r;
    uv = ru*r1;
    vv = rv*r1;
    wv = rw*r1;
    E = rE*r1;
    af = 0.5*(uv*uv+vv*vv+wv*wv);
    p = gam1*(rE -r*af);
    h = E + p*r1;
    c2 = gam*p*r1;
    c = sqrt(c2);
    un = uv*nx + vv*ny + wv*nz;
    
    rlam1 = tanh(1e3*(un+c));
    rlam2 = tanh(1e3*(un-c));
    rlam3 = tanh(1e3*(un));
    
    s1 = 0.5*(rlam1 + rlam2);
    s2 = 0.5*(rlam1 - rlam2);
    
    An = zeros(5,5);
    cc1 = gam1*(s1-rlam3)*af/c2-(s2*un/c);
    cc2 = gam1*s2*af/c-(s1-rlam3)*un;
    An1 = [rlam3+cc1; cc1*uv+cc2*nx; cc1*vv+cc2*ny; cc1*wv+cc2*nz; cc1*h+cc2*un];
    cc1 = -gam1*(s1-rlam3)*uv/c2+(s2*nx/c);
    cc2 = -gam1*s2*uv/c + (s1-rlam3)*nx;
    An2 = [cc1; rlam3+cc1*uv+cc2*nx; cc1*vv+cc2*ny; cc1*wv+cc2*nz; cc1*h+cc2*un];
    cc1 = -gam1*(s1-rlam3)*vv/c2+(s2*ny/c);
    cc2 = -gam1*s2*vv/c+(s1-rlam3)*ny;
    An3 = [cc1; cc1*uv+cc2*nx; rlam3+cc1*vv+cc2*ny; cc1*wv+cc2*nz; cc1*h+cc2*un];
    cc1 = -gam1*(s1-rlam3)*wv/c2+(s2*nz/c);
    cc2 = -gam1*s2*wv/c+(s1-rlam3)*nz;
    An4 = [cc1; cc1*uv+cc2*nx; cc1*vv+cc2*ny; rlam3+cc1*wv+cc2*nz; cc1*h+cc2*un];
    cc1 = gam1*(s1-rlam3)/c2;
    cc2 = gam1*s2/c;
    An5 = [cc1; cc1*uv+cc2*nx; cc1*vv+cc2*ny; cc1*wv+cc2*nz; rlam3+cc1*h+cc2*un];
    
    An = [An1, An2, An3, An4, An5];
    
    % freestream boundary condition
    uinf = sym(mu(3:7)); % freestream flow
    uinf = uinf(:);
    u = u(:);          % state variables
    ui = 0.5*((u+uinf) + An*(u-uinf));  % Riemann solution
    fi = f(:,1)*nx + f(:,2)*ny + f(:,3)*nz + tau*(u-ui); % numerical flux at freestream boundary
    
    %fw = [0; p*nx; p*ny; p*nz; 0]; % inviscid wall boundary condition
    uw = u;
    run = uw(2)*nx+ uw(3)*ny + uw(4)*nz;
    uw(2) = ru - run*nx;
    uw(3) = rv - run*ny;
    uw(3) = rw - run*nz;
    uw = 0.5*((u+uw) + An*(u-uw));  % Riemann solution
    fw = f(:,1)*nx + f(:,2)*ny + f(:,3)*nz + tau*(u-uw); % numerical flux at freestream boundary
    fw(1) = 0;
    fw(end) = 0;
    
    fb = [fw fi];
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    ub = sym(zeros(5,3));
end

function u0 = initu(x, mu, eta)
    u0 = mu(3:7);
end
