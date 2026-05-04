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
m = sym([1.0; 1.0]);
end

function f = flux(u, q, w, v, x, t, mu, eta)
g = mu(1);
Hmin = mu(7);

phi = u(1);
mom = u(2);
H = Hmin + exp(phi);

f = [exp(-phi)*mom; mom*mom/H + 0.5*g*H*H];
end

function s = source(u, q, w, v, x, t, mu, eta)
g = mu(1);
Hmin = mu(7);

phi = u(1);
mom = u(2);
phix = -q(1);
H = Hmin + exp(phi);
[~, dzbdx] = bathymetry(x(1), mu);

s1 = -exp(-phi)*mom*phix;
s2 = -g*H*dzbdx;
s = [s1; s2];
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
uL = leftstate(u, mu);
uR = rightstate(u, x, t, mu);

fbL = rusanovfluxn(u, uL, n(1), mu);
fbR = rusanovfluxn(u, uR, n(1), mu);
fb = [fbL fbR];
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
ub = sym(zeros(2,2));
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
uL = leftstate(u, mu);
uR = rightstate(u, x, t, mu);

fhat = flux(uhat, q, w, v, x, t, mu, eta)*n(1);
fbL = rusanovfluxn(uhat, uL, n(1), mu) - fhat;
fbR = rusanovfluxn(uhat, uR, n(1), mu) - fhat;
fb = [fbL fbR];

% fbL = uL - uhat;
% fbR = uR - uhat;
% fb = [fbL fbR];
end

function u0 = initu(x, mu, eta)
Hmin = mu(7);

%Heps = mu(8);
[zb, ~] = bathymetry(x(1), mu);
% H0 = smoothmax(-zb, Hmin + Heps, Heps);
% phi0 = log(H0 - Hmin);

gamma = 1.0e2;
h = - zb;
H0 = h.*(atan(gamma*(h))/pi + 0.5) - atan(gamma)/pi + 0.5 + Hmin;    
phi0 = log(H0);

u0 = [phi0; sym(0.0)];
end

function uinf = leftstate(u, mu)
Hmin = mu(7);
%Heps = mu(8);
phiL = log(Hmin);
uinf = [phiL; u(2)];
end

function uinf = rightstate(u, x, t, mu)
A = mu(5);
omega = mu(6);
Hmin = mu(7);
%Heps = mu(8);

[zb, ~] = bathymetry(x(1), mu);
%Hbc = smoothmax(A*sin(omega*t) - zb, Hmin + Heps, Heps);
Hbc = A*sin(omega*t) - zb;
Hint = Hmin + exp(u(1));
uint = u(2)/Hint;
phibc = log(Hbc - Hmin);
uinf = [phibc; Hbc*uint];
end

function fh = rusanovfluxn(uL, uR, nx, mu)
fL = flux(uL, [], [], [], [], [], mu, []);
fR = flux(uR, [], [], [], [], [], mu, []);
aL = wavespeed(uL, mu);
aR = wavespeed(uR, mu);
a = smoothmax(aL, aR, 1e-8);

fh = 0.5*(fL + fR)*nx + 0.5*a*(uL(:) - uR(:));
end

function a = wavespeed(u, mu)
g = mu(1);
Hmin = mu(7);

phi = u(1);
mom = u(2);
H = Hmin + exp(phi);
uvel = mom/H;
%a = abs(uvel) + sqrt(g*H);
a = 1 + sqrt(g*H);
end

function y = smoothmax(a, b, eps0)
y = 0.5*(a + b + sqrt((a - b)*(a - b) + eps0*eps0));
end

function [zb, dzbdx] = bathymetry(x, mu)
L = mu(2);
zbeach = mu(3);
zdeep = mu(4);

slope = (zdeep - zbeach)/L;
zb = zbeach + slope*x;
dzbdx = slope;
end
