function pde = pdemodel
pde.mass = @mass;
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.fint = @fint;
pde.ubou = @ubou;
pde.initu = @initu;
% pde.avfield = @avfield;
pde.fbouhdg = @fbouhdg;
end

function m = mass(u, q, w, v, x, t, mu, eta)
m = sym([1.0; 1.0; 1.0; 1.0]); 
end

function f = flux(u, q, w, v, x, t, mu, eta)
   % pde.physicsmu = [gam Re Pr Minf rinf ruinf rvinf rEinf Tref avk avs];

gam = mu(1);
gam1 = gam - 1.0;
Re = mu(2);
Pr = mu(3);
Minf = mu(4);
Tref = mu(10);
muRef = 1/Re;
Tinf = 1/(gam*gam1*Minf^2);
c23 = 2.0/3.0;

% regularization mueters
alpha = 1.0e3;
rmin = 1.0e-2;
pmin = 1.0e-3;

av = v(1);

r = u(1);
ru = u(2);
rv = u(3);
rE = u(4);
rx = q(1);
rux = q(2);
rvx = q(3);
rEx = q(4);
ry = q(5);
ruy = q(6);
rvy = q(7);
rEy = q(8);

% Regularization of rho (cannot be smaller than rmin)
r = rmin + lmax(r-rmin,alpha);
% Density sensor
dr = atan(alpha*(r - rmin))/pi + (alpha*(r - rmin))/(pi*(alpha^2*(r - rmin)^2 + 1)) + 1/2;
%dr=1;
rx = rx*dr;
ry = ry*dr;
r1 = 1/r;
uv = ru*r1;
vv = rv*r1;
E = rE*r1;
q = 0.5*(uv*uv+vv*vv);
p = gam1*(rE-r*q);
% Regularization of pressure p (cannot be smaller than pmin)
p = pmin + lmax(p-pmin,alpha);
% Pressure sensor
dp = atan(alpha*(p - pmin))/pi + (alpha*(p - pmin))/(pi*(alpha^2*(p - pmin)^2 + 1)) + 1/2;
%dp=1;
% Total enthalpy
h = E+p*r1;
% Inviscid fluxes
fi = [ru, ru*uv+p, rv*uv, ru*h, ...
        rv, ru*vv, rv*vv+p, rv*h];
ux = (rux - rx*uv)*r1;
vx = (rvx - rx*vv)*r1;
qx = uv*ux + vv*vx;
px = gam1*(rEx - rx*q - r*qx);
px = px*dp;
Tx = 1/gam1*(px*r - p*rx)*r1^2;
uy = (ruy - ry*uv)*r1;
vy = (rvy - ry*vv)*r1;
qy = uv*uy + vv*vy;
py = gam1*(rEy - ry*q - r*qy);
py = py*dp;
Ty = 1/gam1*(py*r - p*ry)*r1^2;
% Adding Artificial viscosities
T = p/(gam1*r);
Tphys = Tref/Tinf * T;
mu = getViscosity(muRef,Tref,Tphys,1);
mu = mu;
fc = mu*gam/(Pr);
% Viscous fluxes with artificial viscosities
txx = (mu)*c23*(2*ux - vy);
txy = (mu)*(uy + vx);
tyy = (mu)*c23*(2*vy - ux);
fv = [0, txx, txy, uv*txx + vv*txy + (fc)*Tx, ...
      0, txy, tyy, uv*txy + vv*tyy + (fc)*Ty];
fl = [av.*rx, av.*rux, av.*rvx, av.*rEx, av.*ry, av.*ruy, av.*rvy, av.*rEy];
f = fi+fv +fl;

    f = reshape(f,[4,2]);        
end

% function f = avfield(u, q, w, v, x, t, mu, eta)
%     f = getavfield2d(u,q,v,mu);
% end

function s = source(u, q, w, v, x, t, mu, eta)
    s = [sym(0.0); sym(0.0); sym(0.0); sym(0.0)];
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
   
    fb = [sym(0.0); sym(0.0); sym(0.0); sym(0.0)];
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    ub = [sym(0.0); sym(0.0); sym(0.0); sym(0.0)];
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)


    gam = mu(1);
    gam1 = gam - 1.0;
    Tinf = mu(9);
    Tref = mu(10);
    Twall = mu(11);
    TisoW = Twall/Tref * Tinf;    
    uinf = sym(mu(5:8)); % freestream flow
    uinf = uinf(:);

    f_out = u - uhat;
    f_in = uinf - uhat;

    % iso-thermal wall boundary condition    
    f1 = 0*u;
    f1(1) = u(1) - uhat(1); % extrapolate density
    f1(2) = 0.0  - uhat(2); % zero velocity
    f1(3) = 0.0  - uhat(3); % zero velocity           
    f1(4) = -uhat(4) + uhat(1)*TisoW;
    
    % adiabatic wall boundary condition    
    f = flux(uhat, q, w, v, x, t, mu, eta);
    f2 = 0*u;
    f2(1) = u(1) - uhat(1); % extrapolate density
    f2(2) = 0.0  - uhat(2); % zero velocity
    f2(3) = 0.0  - uhat(3); % zero velocity               
    f2(4) = f(4,1)*n(1) + f(4,2)*n(2) + tau*(u(4)-uhat(4)); % zero heat flux
    
    % iso-thermal wall boundary condition    
    f3 = 0*u;
    f3(1) = u(1) - uhat(1); % extrapolate density
    f3(2) = 0.0  - uhat(2); % zero velocity
    f3(3) = 0.0  - uhat(3); % zero velocity           
    f3(4) = uhat(4)/uhat(1);
    
    % supersonic inflow, supersonic outflow, isothermal, adiabatic, isothermal
    fb = [f_in f_out f1 f2 f3];
end

function u0 = initu(x, mu, eta)
    u0 = sym(mu(5:8)); % freestream flow   
end

function fi = fint(u, q, w, v, x, t, mu, eta, uhat, n, tau)
f = flux(uhat, q, w, v, x, t, mu, eta);
fi = f(4,1)*n(1) + f(4,2)*n(2) + tau*(u(4)-uhat(4)); 
end


