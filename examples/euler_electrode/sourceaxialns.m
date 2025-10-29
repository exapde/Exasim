function f = sourceaxialns(u, q, w, v, x, t, mu, eta)
   % pde.physicsmu = [gam Re Pr Minf rinf ruinf rvinf rEinf Tref avk avs];

% regularization mueters
alpha = 4.0e3;
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
% Density sensor
dr = atan(alpha*(r - rmin))/pi + (alpha*(r - rmin))/(pi*(alpha^2*(r - rmin)^2 + 1)) + 1/2;
r1 = 1/r;
uv = ru*r1;
vv = rv*r1;
E = rE*r1;

gam = 1.4;        % Thermodynamic quantities are determined consistently from the equation of state and used to couple the energy to the momentum equation in the following manner: energy->temperature->pressure
e = E-0.5*(uv*uv+vv*vv);
T = e*(gam-1);
H = E + T;

% Everything above is same as the flux term, with the exception that the pressure has been removed.


% Inviscid fluxes
fi = [rv, ru*vv, rv*vv, rv*H];

f = -fi/x(2);  % For axisymmetric, radial coordinate is x(2)

f = reshape(f,[4,1]);        

end

