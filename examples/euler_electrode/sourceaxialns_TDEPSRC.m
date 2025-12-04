function f = sourceaxialns_TDEPSRC(u, q, w, v, x, t, mu, eta)
   % pde.physicsmu = [gam Re Pr Minf rinf ruinf rvinf rEinf Tref avk avs];

   disp("SOURCE TERM TDEP")

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

% Time dependent source term
tfinal = 50;  % Final time for distribution to be "over" in nondimensional simulation time
t_target = 100e-6;   % Final (target) time to plot over. For the heat release from a gas kernel, this should be on the order of tens of microseconds
x_scale = tfinal/t_target;       % Assume a linear scaling in x for now, but you could generalize this to the nonlinear case
total_energy = 300e-6;
k=1.5;
Lref = 1e-4;
Pref = 101325;
rho_ref = 1.225;
t_ref = Lref/(sqrt(Pref/rho_ref));   % = 3.4770e-07
spark_gap_length = 1.5e-3;  % m
spark_kernel_radius = 100e-6;  % m
spark_kernel_volume = pi*spark_kernel_radius^2*spark_gap_length;  % m
Eref = (Pref*sqrt(Pref/rho_ref))/Lref;

phys_t = t*t_ref;  % t scale
energy_rate_dimensional = total_energy*x_scale*k*exp(-k*x_scale*phys_t);    % J/s - PHYSICAL VALUE

% divide by geometry to get the energy DENSITY deposition rate, J/m^3-s
energy_density_rate_dimensional = energy_rate_dimensional/spark_kernel_volume;       % J/(m^3-s) - PHYSICAL VALUE - same units as the energy conservative variable

energy_density_rate_temporal_NON_dimensional = energy_density_rate_dimensional / Eref;
energy_density_rate_temporal_spatial_NON_dimensional = energy_density_rate_temporal_NON_dimensional/(1+exp(20*(x(2)-1)));    % Scale by sigmoid to get the spatial profile

energy_term = rv*H - energy_density_rate_temporal_spatial_NON_dimensional;

fi = [rv, ru*vv, rv*vv, energy_term];

av = v(1);
fl = [av*ry, av*ruy, av*rvy, av*rEy];

f = fi + fl;

f = -f/(0.1+x(2));  % For axisymmetric, radial coordinate is x(2)

f = reshape(f,[4,1]);        

end
