function pde = pdemodel
pde.mass = @mass;
pde.flux = @flux;
pde.source = @source;
pde.sourcew = @sourcew;
pde.fbou = @fbou;
pde.ubou = @ubou;
pde.initu = @initu;
% pde.avfield = @avfield;
pde.fbouhdg = @fbouhdg;
pde.fhathdg = @fhathdg;
pde.visscalars = @visscalars;
pde.visvectors = @visvectors;
end

function m = mass(u, q, w, v, x, t, mu, eta)
m = sym([1.0; 1.0; 1.0; 1.0]);
end

function f = flux(u, q, w, v, x, t, mu, eta)

f = fluxaxialns_ORIG(u, q, w, v, x, t, mu, eta);
% f = fluxaxialns_SOURCEW(u, q, w, v, x, t, mu, eta);

end

% function f = avfield(u, q, w, v, x, t, mu, eta)
%     f = getavfield2d(u,q,v,mu);
% end

function s = visscalars(u, q, w, v, x, t, mu, eta)
s(1) = u(1);
s(2) = u(4);






% % Sponge layer has two paramters: length and strength. https://web.stanford.edu/group/ctr/ResBriefs/2010/11_mani.pdf
% r0 = 50;   % nondimensional sponge start
% beta = 1e-2;
% sigma = -beta*(x(2)-r0)^2*switchSigmoid(x(2), r0, 1);

% % s(1) = x(2);
% % s(2) = switchSigmoid(x(2), r0, 1);
% % s(3) = (x(2)-r0)^2;
% % s(4) = s(2)*s(3);
% % s(5) = -beta*s(4);

% % % Farfield nondimensional values for the conservative variables
% rho_0 = 1.225;
% ru0 = 0;
% rv0 = 0;
% rE0 = 2.5;

% s(3) = sigma*(u(1)-rho_0);
% s(4) = sigma*(u(2)-ru0);
% s(5) = sigma*(u(3)-rv0);
% s(6) = sigma*(u(4)-rE0);


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

s(3) = t;    % Simulation time
s(4) = phys_t;     % Physical time
s(5) = energy_density_rate_temporal_spatial_NON_dimensional;





% s(3) = q(1);
% s(4) = q(5);

% s(5) = q(2);
% s(6) = q(6);

% s(7) = q(3);
% s(8) = q(7);

% s(9) = q(4);
% s(10) = q(8);
end

function s = visvectors(u, q, w, v, x, t, mu, eta)
s(1) = u(2);
s(2) = u(3);



% s(3) = q(1);
% s(4) = q(5);

% s(5) = q(2);
% s(6) = q(6);

% s(7) = q(3);
% s(8) = q(7);

% s(9) = q(4);
% s(10) = q(8);
end

function s = source(u, q, w, v, x, t, mu, eta)

% s = sourceaxialns_SPONGE(u, q, w, v, x, t, mu, eta);
s = sourceaxialns_TDEPSRC(u, q, w, v, x, t, mu, eta);

end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)

fb = [sym(0.0); sym(0.0); sym(0.0); sym(0.0)];
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
ub = [sym(0.0); sym(0.0); sym(0.0); sym(0.0)];
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)

fb = fbouhdgaxialns(u, q, w, v, x, t, mu, eta, uhat, n, tau);

end

function fh = fhathdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)

f = flux(uhat, q, w, v, x, t, mu, eta);

r = uhat(1);
ru = uhat(2);
rv = uhat(3);
nx = n(1);
ny = n(2);
run = ru*nx + rv*ny;
r1 = 1/r;
un = run*r1;

s = (u - uhat);
for i = 1:4
    s(i) = r*tau(i)*(u(i) - uhat(i));
end
% s(2) = sqrt(un*un + 1e-10)*(u(2) - uhat(2));
% s(3) = sqrt(un*un + 1e-10)*(u(3) - uhat(3));

fh = f(:,1)*n(1) + f(:,2)*n(2) + s;

end

function u0 = initu(x, mu, eta)
u0 = [sym(0.0); sym(0.0); sym(0.0); sym(0.0)];
end


