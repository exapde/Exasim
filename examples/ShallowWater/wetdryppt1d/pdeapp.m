% Add Exasim to Matlab search path.
cdir = pwd();
ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

% Initialize PDE and mesh structures.
[pde, mesh] = initializeexasim();

% Positivity-preserving transformed shallow-water model:
% H(phi) = Hmin + exp(phi), with ModelD used for auxiliary gradients.
pde.model = "ModelD";
pde.modelfile = "pdemodel";

% Use standard DG/LDG for the transformed system.
pde.platform = "cpu";
pde.mpiprocs = 1;
pde.hybrid = 1;
pde.Cxxpreprocessing = 0;

% Physical parameters.
g = 9.81;
L = 10000.0;
zbeach = 3.0;
zdeep = -7.0;
A_tide = 0.4;
T_tide = 1800.0;
omega = 2.0*pi/T_tide;
Hmin = 1e-3;
Heps = 1e-8;
pde.physicsparam = [g, L, zbeach, zdeep, A_tide, omega, Hmin];

% Discretization and time integration.
N = 400;
dt = 1.0;
nsteps = 100;

pde.porder = 1;
pde.torder = 2;
pde.nstage = 2;
pde.tau = 4;
pde.dt = dt*ones(1, nsteps);
pde.saveSolFreq = 1;
pde.soltime = pde.saveSolFreq:pde.saveSolFreq:nsteps;
pde.visdt = pde.saveSolFreq*dt;

pde.GMRESrestart = 30;
pde.linearsolvertol = 1e-8;
pde.linearsolveriter = 30;
pde.NLtol = 1e-8;
pde.NLiter = 3;

% Uniform 1D mesh on [0, L].
[mesh.p, mesh.t] = linemesh(N);
mesh.p = L*mesh.p;
mesh.boundaryexpr = {@(p) abs(p(1,:)) < 1e-8, @(p) abs(p(1,:) - L) < 1e-8};
mesh.boundarycondition = [1; 2];

% Generate and run the Exasim model.
[sol, pde, mesh, master, dmd] = exasim(pde, mesh);

% Postprocess in physical variables.
dgnodes = createdgnodes(mesh.p, mesh.t, mesh.f, mesh.curvedboundary, ...
    mesh.curvedboundaryexpr, pde.porder);
x = squeeze(dgnodes(:,1,:));
phi = squeeze(sol(:,1,:,end));
mom = squeeze(sol(:,2,:,end));

zb = localbathymetry(x, pde.physicsparam);
H = Hmin + exp(phi);
eta = H + zb;
u = mom./H;

for i = 1:size(sol,4)
    phi = squeeze(sol(:,1,:,i));
    mom = squeeze(sol(:,2,:,i));
    H = Hmin + exp(phi);
    eta = H + zb;
    figure(2); clf; plot(x(:), H(:), "LineWidth", 1.2); axis([0 max(x(:)) 0 8]);
    xlabel("x [m]");
    ylabel("H [m]");
    title("Transformed depth evolution");
    grid on;
end

xplot = x(:);
phiplot = phi(:);
Hplot = H(:);
etaplot = eta(:);
uplot = u(:);
zbplot = zb(:);
[xplot, ind] = sort(xplot);
phiplot = phiplot(ind);
Hplot = Hplot(ind);
etaplot = etaplot(ind);
uplot = uplot(ind);
zbplot = zbplot(ind);

figure(1); clf;
subplot(4,1,1);
plot(xplot, phiplot, "LineWidth", 1.5);
xlabel("x [m]");
ylabel("\phi");
title("Transformed variable");
grid on;

subplot(4,1,2);
plot(xplot, Hplot, "LineWidth", 1.5);
xlabel("x [m]");
ylabel("H [m]");
title("Water depth");
grid on;

subplot(4,1,3);
plot(xplot, etaplot, "b-", "LineWidth", 1.5); hold on;
plot(xplot, zbplot, "k--", "LineWidth", 1.0);
xlabel("x [m]");
ylabel("elevation [m]");
legend("free surface", "bed", "Location", "best");
title("Free surface and bathymetry");
grid on;

subplot(4,1,4);
plot(xplot, uplot, "LineWidth", 1.5);
xlabel("x [m]");
ylabel("u [m/s]");
title("Depth-averaged velocity");
grid on;

disp("Done!");

