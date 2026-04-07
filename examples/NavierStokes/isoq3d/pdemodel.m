function pde = pdemodel
pde.mass = @mass;
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.fbouhdg = @fbouhdg;
pde.ubou = @ubou;
pde.initu = @initu;
end

function m = mass(u, q, w, v, x, t, mu, eta)
m = sym([1.0; 1.0; 1.0; 1.0; 1.0]); 
end

function f = flux(u, q, w, v, x, t, mu, eta)
    gam = mu(1);    
    gam1 = gam - 1.0;
    Re = mu(2);
    Pr = mu(3);
    Minf = mu(4);

    Tinf = mu(10);
    Tref = mu(11);
    %Twall = mu(12);
    muRef = 1/Re;
    Re1 = 1/Re;
    M2 = Minf^2;
    c23 = 2.0/3.0;

    % regularization params
    alpha = 1.0e3;
    rmin = 1.0e-2;
    pmin = 1.0e-3;

    av = v(1);

    %fc = 1/(gam1*M2*Re*Pr);
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

    % Regularization of rho (cannot be smaller than rmin)
    r = rmin + lmax(r-rmin,alpha);
    % Density sensor
    dr = atan(alpha*(r - rmin))/pi + (alpha*(r - rmin))/(pi*(alpha^2*(r - rmin)^2 + 1)) + 1/2;
    rx = rx*dr;
    ry = ry*dr;
    rz = rz*dr;

    r1 = 1/r;
    uv = ru*r1;
    vv = rv*r1;
    wv = rw*r1;
    E = rE*r1;
    ke = 0.5*(uv*uv+vv*vv+wv*wv);
    p = gam1*(rE-r*ke);

    % Regularization of pressure p (cannot be smaller than pmin)
    p = pmin + lmax(p-pmin,alpha);
    % Pressure sensor
    dp = atan(alpha*(p - pmin))/pi + (alpha*(p - pmin))/(pi*(alpha^2*(p - pmin)^2 + 1)) + 1/2;

    % Total enthalpy
    h = E+p*r1;
    % Inviscid fluxes
    fi = [ru, ru*uv+p, rv*uv, rw*uv, ru*h, ...
          rv, ru*vv, rv*vv+p, rw*vv, rv*h, ...
          rw, ru*wv, rv*wv, rw*wv+p, rw*h];


    %x-component
    ux = (rux - rx*uv)*r1;
    vx = (rvx - rx*vv)*r1;
    wx = (rwx - rx*wv)*r1;
    qx = uv*ux + vv*vx + wv*wx;
    px = gam1*(rEx - rx*ke - r*qx);
    px = px*dp;
    Tx = 1/gam1*(px*r - p*rx)*r1^2;

    %y-component
    uy = (ruy - ry*uv)*r1;
    vy = (rvy - ry*vv)*r1;
    wy = (rwy - ry*wv)*r1;
    qy = uv*uy + vv*vy + wv*wy;
    py = gam1*(rEy - ry*ke - r*qy);
    py = py*dp;
    Ty = 1/gam1*(py*r - p*ry)*r1^2;

    %z-component
    uz = (ruz - rz*uv)*r1;
    vz = (rvz - rz*vv)*r1;
    wz = (rwz - rz*wv)*r1;
    qz = uv*uz + vv*vz + wv*wz;
    pz = gam1*(rEz - rz*ke - r*qz);
    pz = pz*dp;
    Tz = 1/gam1*(pz*r - p*rz)*r1^2;

    % Adding Artificial viscosities
    T = p/(gam1*r);
    Tphys = Tref/Tinf * T;
    mu = getViscosity(muRef,Tref,Tphys,1);
    fc = mu*gam/(Pr);

    % Viscous fluxes with artificial viscosities
    txx = (mu)*c23*(2*ux - vy - wz);
    txy = (mu)*(uy + vx);
    txz = (mu)*(uz + wx);
    tyy = (mu)*c23*(2*vy - ux - wz);
    tyz = (mu)*(vz + wy);
    tzz = (mu)*c23*(2*wz - ux - vy);
    fv = [0, txx, txy, txz, uv*txx + vv*txy + wv*txz + (fc)*Tx,...
          0, txy, tyy, tyz, uv*txy + vv*tyy + wv*tyz + (fc)*Ty,...
          0, txz, tyz, tzz, uv*txz + vv*tyz + wv*tzz + (fc)*Tz];
    %artifical viscosities
    fl = [av.*rx, av.*rux, av.*rvx, av.*rwx, av.*rEx,...
          av.*ry, av.*ruy, av.*rvy, av.*rwy, av.*rEy,...
          av.*rz, av.*ruz, av.*rvz, av.*rwz, av.*rEz];

    f = fi+fv+fl;
    f = reshape(f,[5,3]);    
end

function s = source(u, q, w, v, x, t, mu, eta)
s = [sym(0.0); sym(0.0); sym(0.0); sym(0.0); sym(0.0)];
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
f = flux(u, q, w, v, x, t, mu, eta);
fi = f(:,1)*n(1) + f(:,2)*n(2) + f(:,3)*n(3) + tau*(u-uhat);

% thermal Wall
ftw = fi;
ftw(1) = 0.0;

% inflow, outflow, thermal, periodic
fb = [fi fi ftw fi]; 
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)    
% Isothermal Wall
Tinf = mu(10);
Tref = mu(11);
Twall = mu(12);
TisoW = Twall/Tref * Tinf;    
utw = u(:);
utw(2:4) = 0;
utw(5) = u(1)*TisoW;
    
% inflow, outflow, thermal wall, periodic
ub = [u utw  0*u]; 
end

function u0 = initu(x, mu, eta)
    u0 = sym(mu(5:9)); % freestream flow   
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    gam = mu(1);
    gam1 = gam - 1.0;
    Tinf = mu(10);
    Tref = mu(11);
    Twall = mu(12);
    TisoW = Twall/Tref * Tinf;    
    uinf = sym(mu(5:9)); % freestream flow
    uinf = uinf(:);

    f_out = u - uhat;
    f_in = uinf - uhat;

    % wall boundary condition    
    f1 = 0*u;
    f1(1) = u(1) - uhat(1); % extrapolate density
    f1(2) = 0.0  - uhat(2); % zero velocity %x
    f1(3) = 0.0  - uhat(3); % zero velocity %y
    f1(4) = 0.0  - uhat(4); % zero velocity %z
    f1(5) = -uhat(5) +uhat(1)*TisoW;
    f_wall = f1;

    % slip wall condition                
    ru = u(2);
    rv = u(3);
    rw = u(4);
    nx = n(1);
    ny = n(2);
    nz = n(3);
    run = ru*nx + rv*ny + rw*nz;        
    uinf = u;
    uinf(2) = uinf(2) - nx.*run;
    uinf(3) = uinf(3) - ny.*run;
    uinf(4) = uinf(4) - nz.*run;
    f_slip = (uinf - uhat);    
    
    % zero gradient condition
    q = q(:);
    f_grad = q(1:5)*n(1) + q(6:10)*n(2) + q(11:15)*n(3) + tau*(u(:) - uhat(:)); 

    %periodic 
    f = flux(uhat, q, w, v, x, t, mu, eta);
    f_periodic = f(:,1)*n(1) + f(:,2)*n(2) + f(:,3)*n(3) + tau*(u-uhat);        
    
    fb = [f_in f_out f_wall f_grad f_periodic f_slip];
end
