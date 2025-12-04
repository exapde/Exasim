function f = fluxaxialns_ORIG(u, q, w, v, x, t, mu, eta)
   % pde.physicsmu = [gam Re Pr Minf rinf ruinf rvinf rEinf Tref avk avs];
      disp("FLUX ORIG")

% regularization mueters
alpha = 1.0e2;
rmin = 5.0e-2;
pmin = 2.0e-3;

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
% % Density sensor
% dr = atan(alpha*(r - rmin))/pi + (alpha*(r - rmin))/(pi*(alpha^2*(r - rmin)^2 + 1)) + 1/2;
r1 = 1/r;
uv = ru*r1;
vv = rv*r1;
E = rE*r1;

gam = 1.4;        % Thermodynamic quantities are determined consistently from the equation of state and used to couple the energy to the momentum equation in the following manner: energy->temperature->pressure
e = E-0.5*(uv*uv+vv*vv);
T = e*(gam-1);
p = r*T;
H = E + T;

%p = pmin + lmax(p-pmin,alpha);

% Inviscid fluxes only -- Euler
fi = [ru, ru*uv+p, rv*uv, ru*H, ...
      rv, ru*vv, rv*vv+p, rv*H];

av = v(1);
fl = [av*rx, av*rux, av*rvx, av*rEx, av*ry, av*ruy, av*rvy, av*rEy];    

f = fi + fl;

f = reshape(f,[4,2]);        

end
