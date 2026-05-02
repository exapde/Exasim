function pde = pdemodel
pde.mass = @mass;
pde.flux = @flux;
pde.source = @source;
pde.avfield = @avfield;
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
hfloor = mu(7);
H = positivepart(u(1));
Hsafe = H + hfloor;
qmom = u(2);

% Hx = q(1);
% Hux = q(2);
% % Hm = 1e-3 - H;
% gamma = 1e3;
% Hsafe = u(1).*(atan(gamma*(u(1)))/pi + 0.5) - atan(gamma)/pi + 0.5 + 1e-4;    
% vel = qmom/Hsafe;
% av = vel.*(atan(gamma*(vel))/pi + 0.5) - atan(gamma)/pi + 0.5;    
% fl = 0*[av.*Hx; av.*Hux];

av = v(1);
Hx = q(1);
Hux = q(2);
fl = [30*av.*Hx; 30*av.*Hux];

f = [qmom; qmom*qmom/Hsafe + 0.5*g*H(1)*H(1)] + fl;
end

function s = source(u, q, w, v, x, t, mu, eta)
g = mu(1);
%hfloor = mu(7);
H = positivepart(u(1));
[~, dzbdx] = bathymetry(x(1), mu);

% Smoothly damp momentum in nearly dry cells while leaving wet cells almost unchanged.
Htol = 1e-3;
deltaH = 1e-3;
alpha0 = 0.1;
chi = 0.5*(1.0 - tanh((H - Htol)/deltaH));
sdamp = -alpha0*chi*u(2);

s = [sym(0.0); -g*H*dzbdx + sdamp];
end

function av = avfield(u, q, w, v, x, t, mu, eta)
  
  Hx = q(1);
  qx = q(2);
  H = u(1);
  q = u(2);
  hfloor = mu(7);
  L = mu(2);
  dist = (L-x)/L;  

  gamma = 1e3;
  Hp = H.*(atan(gamma*(H))/pi + 0.5) - atan(gamma)/pi + 0.5;    
  dHpdH = (atan(gamma*H)/pi + 0.5) + H.*(gamma./(pi*(1 + gamma^2*H.^2)));
  Hpx = dHpdH .* Hx;
  den = (Hp + hfloor);
  ux = (qx.*den - q.*Hpx) ./ (den.^2);
  
  s = ux;
  s = s/0.001;
  s = limiting(s,0.0,0.9,gamma,0.0);
  s = s/0.9;
  S0 = 0.1; 
  av = (s-S0).*(atan(gamma*(s-S0))/pi + 0.5) - atan(gamma)/pi + 0.5; 
  av = av*tanh(10*dist);
end


function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
f = flux(u, q, w, v, x, t, mu, eta);
fb = f*n(1) + tau.*(u(:) - uhat(:));
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
ub = sym(zeros(2,2));
end

function u0 = initu(x, mu, eta)
[zb, ~] = bathymetry(x(1), mu);
H0 = positivepart(-zb);
u0 = [H0; sym(0.0)];
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
u = u(:);
uhat = uhat(:);

uinfL = leftstate(u, mu);
uinfR = rightstate(u, x, t, mu);
An = swematrix(uhat, n(1), mu);

fbL = 0.5*((u + uinfL) + An*(u - uinfL)) - uhat;
fbR = 0.5*((u + uinfR) + An*(u - uinfR)) - uhat;
fb = [fbL fbR];
end

function uinf = leftstate(u, mu)
uinf = [0; 0];
end

function uinf = rightstate(u, x, t, mu)
A = mu(5);
omega = mu(6);

[zb, ~] = bathymetry(x(1), mu);
Hbc = A*sin(omega*t) - zb;
uint = u(2)/u(1);
uinf = [Hbc; Hbc*uint];
end

function An = swematrix(uhat, nx, mu)
g = mu(1);
Hhat = positivepart(uhat(1));
Hsafe = Hhat;
qhat = uhat(2);
uvel = qhat/Hsafe;
c = sqrt(g*Hsafe);

lam1 = nx*(uvel - c);
lam2 = nx*(uvel + c);
beta = 1e2;
Lambda = [tanh(beta*lam1), sym(0.0); sym(0.0), tanh(beta*lam2)];

K = [sym(1.0), sym(1.0); uvel - c, uvel + c];
Kinv = (1/(2*c))*[uvel + c, -sym(1.0); c - uvel, sym(1.0)];

An = simplify(K*Lambda*Kinv);
end

function H = positivepart(h)
gamma = 1e6;
H = h.*(atan(gamma*(h))/pi + 0.5) - atan(gamma)/pi + 0.5;    
end

function [zb, dzbdx] = bathymetry(x, mu)
L = mu(2);
zbeach = mu(3);
zdeep = mu(4);

slope = (zdeep - zbeach)/L;
zb = zbeach + slope*x;
dzbdx = slope;
end
