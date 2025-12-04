function f = sourceaxialns_SPONGE(u, q, w, v, x, t, mu, eta)
   % pde.physicsmu = [gam Re Pr Minf rinf ruinf rvinf rEinf Tref avk avs];
disp("SOURCE TERM SPONGE")

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
H = E + T;

% Everything above is same as the flux term, with the exception that the pressure has been removed.

% Sponge layer has two paramters: length and strength. https://web.stanford.edu/group/ctr/ResBriefs/2010/11_mani.pdf
r0 = 50;   % nondimensional sponge start
beta = 1e-1;
radius = sqrt(x(1)^2 + x(2)^2);
sigma = beta*(radius-r0)^2*switchSigmoid(radius, r0, 1);

% Farfield nondimensional values for the conservative variables
rho_0 = 1.225;
ru0 = 0;
rv0 = 0;
rE0 = 2.5;

s1 = sigma*(u(1)-rho_0);
s2 = sigma*(u(2)-ru0);
s3 = sigma*(u(3)-rv0);
s4 = sigma*(u(4)-rE0);

energy_term = rv*H + s4;
fi = [rv + s1, ru*vv + s2, rv*vv + s3, energy_term];

av = v(1);
fl = [av*ry, av*ruy, av*rvy, av*rEy];

f = fi + fl;

f = -f/(0.1+x(2));  % For axisymmetric, radial coordinate is x(2)

f = reshape(f,[4,1]);        

end

