function pde = pdemodel
pde.mass = @mass;
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.ubou = @ubou;
pde.initu = @initu;
% pde.avfield = @avfield;
pde.fbouhdg = @fbouhdg;
end

function m = mass(u, q, w, v, x, t, mu, eta)
m = sym([1.0; 1.0; 1.0; 1.0; 1.0]); 
end

function f = flux(u, q, w, v, x, t, mu, eta)
% pde.physicsparam = [gam Re Pr Minf rinf ruinf rvinf rEinf Tref avk avs];

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
rN   = u(5);

rx = q(1);
rux = q(2);
rvx = q(3);
rEx = q(4);
rNx = q(5);

ry = q(6);
ruy = q(7);
rvy = q(8);
rEy = q(9);
rNy = q(10);

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
N = rN*r1;
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
fi = [ru, ru*uv+p, rv*uv, ru*h, rN*uv ...
      rv, ru*vv, rv*vv+p, rv*h, rN*vv];

ux = (rux - rx*uv)*r1;
vx = (rvx - rx*vv)*r1;
qx = uv*ux + vv*vx;
px = gam1*(rEx - rx*q - r*qx);
px = px*dp;
Tx = 1/gam1*(px*r - p*rx)*r1^2;
Nx  = (rNx - rx*N)*r1;

uy = (ruy - ry*uv)*r1;
vy = (rvy - ry*vv)*r1;
qy = uv*uy + vv*vy;
py = gam1*(rEy - ry*q - r*qy);
py = py*dp;
Ty = 1/gam1*(py*r - p*ry)*r1^2;
Ny = (rNy - ry*N)*r1;

% molecular viscosity
T = p/(gam1*r);
Tphys = Tref/Tinf * T;
muM = getViscosity(muRef,Tref,Tphys,1);

% turbulent eddy viscocity
cv1  = 7.1;
b = 100;
c = 1/2-atan(b)/pi;
chi = rN*Re./r;
shi = chi.*(atan(b*chi)/pi+1/2)+c; 
fv1 = shi.^3./(shi.^3+cv1^3);
muN = muRef*r.*shi.*fv1;  

% total viscosity
muT = muM + muN;
fc = muT*gam/(Pr);

muS = muRef*shi*r + muRef*r;
fs  = 1.5*muS;

% Viscous fluxes 
txx = (muT)*c23*(2*ux - vy);
txy = (muT)*(uy + vx);
tyy = (muT)*c23*(2*vy - ux);

fv = [0, txx, txy, uv*txx + vv*txy + (fc)*Tx, fs*Nx ...
      0, txy, tyy, uv*txy + vv*tyy + (fc)*Ty, fs*Ny];

% artificial viscosity fluxes
fl = [av*rx, av*rux, av*rvx, av*rEx, av*rNx, av*ry, av*ruy, av*rvy, av*rEy, av*rNy];
f = fi+fv +fl;

f = reshape(f,[5,2]);        
end

function s = source(u, q, w, v, x, t, mu, eta)

Re = mu(2);
b = 100.0;
sigma = 2.0/3.0;
cv1   = 7.1;
cb1   = 0.1355;
cb2   = 0.622;
kappa = 0.41;
% cw1   = 3.2391;
cw1   = cb1/kappa^2 + (1+cb2)/sigma;
cw2   = 0.3;
cw3   = 2.0;
rlim  = 10;

r = u(1);
ru = u(2);
rv = u(3);
rN   = u(5);

rx = q(1);
rvx = q(3);
rNx = q(5);

ry = q(6);
ruy = q(7);
rNy = q(10);

u  = ru/r;
v  = rv/r;
N  = rN/r;              % \tilde{nu}, eddy viscosity
nu = 1/Re;              % molecular viscosity
vx = (rvx - rx*v)/r;    % will actually be equal to -dv/dx in the code.
uy = (ruy - ry*u)/r;    % will actually be equal to -du/dy in the code.
Nx = (rNx - rx*N)/r;    % will actually be equal to -dN/dx in the code.
Ny = (rNy - ry*N)/r;    % will actually be equal to -dN/dy in the code.

chi = N/nu;
c   = 1/2-atan(b)/pi;
psi = chi*(atan(b*chi)/pi+1/2)+c;   % atan-regularized version of chi

Ds  = (cb2/sigma)*(rNx*rNx+rNy*rNy)/r - (1/sigma)*(nu+N)*(rx*Nx+ry*Ny);    % Diffusion/propagation term (LJ3-p76)

S      = ((uy-vx)*(uy-vx) + 1.0e-20)^(0.5);   % NOTE: (uy-vx) = -omega, not +omega. (Because HDG gradients are stored as -gradients.) But this is ok here.
fv1    = psi^3/(psi^3+cv1^3);
fv2    = 1.0-psi/(1.0+psi*fv1);
Sbar   = (psi*nu)/(kappa^2*d^2)*fv2;
Stilde = 0.1*S + (Sbar+0.9*S)*(atan(b*(Sbar/S+0.9))/pi+1/2)+c*S;  % atan-regularized version of Stilde
Sp     = cb1*Stilde*(r*psi*nu);                  % Production term (note the r, should be there in compressible formulation -- Allmaras 2012)

rc = (psi*nu)/(Stilde*kappa^2*d^2);
rd = rlim-(rlim-rc)*(atan(b*(rlim-rc))/pi+1/2)-c;       % atan-regularized version of rd
gd = rd + cw2*(rd^(6.0) - rd);
fw = gd*((1+cw3^6)/(gd^6 + cw3^6))^(1.0/6.0);
Sd = cw1*fw*r*((psi*nu)/d)^2;                  % Destruction term (note the r; consistent with compressible formulation, Allmaras 2012)

sa = Ds + Sp - Sd;

s = [sym(0.0); sym(0.0); sym(0.0); sym(0.0); sa];    
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    fb = 0*uhat; 
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    ub = 0*uhat; 
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)

gam = mu(1);
Minf = mu(4);
gam1 = gam - 1.0;
r = uhat(1);
ru = uhat(2);
rv = uhat(3);
rE = uhat(4);
nx = n(1);
ny = n(2);

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

K = [ 1 , 1 , 0 , 1 ;...
      un-a , un , 0 , un+a ;...
      ut , ut , 1 , ut ;...
      h - un*a , (1/2)*(un^2 + ut^2) , ut , h+un*a ];
Kinv = (gam1/(2*a^2))*[ h + (a/gam1)*(un-a) , -(un+a/gam1) , -ut , 1 ;...
                        -2*h + (4/gam1)*a^2 , 2*un , 2*ut , -2 ;...
                        -2*(ut*a^2)/gam1 , 0 , 2*(a^2)/gam1 , 0 ;...
                        h - a*(un+a)/gam1 , -un+a/gam1 , -ut , 1 ];
T = [ 1 , 0 , 0 , 0;...
      0 , nx , ny , 0;...
      0 , -ny , nx , 0;...
      0 , 0 , 0 , 1];
Tinv = [ 1 , 0 , 0 , 0;...
         0 , nx ,-ny , 0;...
         0 , ny , nx , 0;...
         0 , 0 , 0 , 1];
Lambda = [ tanh(1e2*(un-a)) , 0 , 0 , 0 ;...
                 0 , tanh(1e2*(un)) , 0 , 0 ;...
                 0 , 0 , tanh(1e2*(un)) , 0 ;...
                 0 , 0 , 0 , tanh(1e2*(un+a)) ];

L = simplify(Tinv * K);
R = simplify(Kinv * T);
An = simplify(L * Lambda * R);

% inflow boundary condition
uinf = mu(5:8);      % freestream flow variables   
u = u(1:4);          % state variables 
fns = 0.5*((u(:)+uinf(:)) + An*(u(:)-uinf(:))) - uhat(:);     
fb1 = [fns; mu(9)-uhat(5)];

% wall boundary condition    
fb2 = 0*fb1;
fb2(1) = u(1) - uhat(1); % extrapolate density
fb2(2) = 0.0  - uhat(2); % zero velocity
fb2(3) = 0.0  - uhat(3); % zero velocity           
f = flux(uhat, q, w, v, x, t, mu, eta);
fb2(4) = f(4,1)*n(1) + f(4,2)*n(2) + tau*(u(4)-uhat(4)); % zero heat flux
fb2(5) = 0.0  - uhat(5); % zero eddy viscosity

% slip wall condition                
ru = u(2);
rv = u(3);
nx = n(1);
ny = n(2);   
run = ru*nx + rv*ny;        
uinf = u;
uinf(2) = uinf(2) - nx.*run;
uinf(3) = uinf(3) - ny.*run;
fb3 = tau*(uinf - uhat);    

% outflow boundary condition  
pinf = 1/(gam*Minf^2); 
uinf = [u(1), u(2), u(3), pinf/(gam-1) + 0.5*u(2)*u(2)/u(1) + 0.5*u(3)*u(3)/u(1)];
u = u(1:4);          % state variables 
fns = 0.5*((u(:)+uinf(:)) + An*(u(:)-uinf(:))) - uhat(:);     
fb4 = [fns; u(5)-uhat(5)];

fb = [fb1 fb2 fb3 fb4];    
end

function u0 = initu(x, mu, eta)
    u0 = sym(mu(5:8)); % freestream flow   
end




