% Add Exasim to Matlab search path.
cdir = pwd();
ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

% Initialize PDE and mesh structures.
[pde, mesh] = initializeexasim();

% Governing equations: 1D shallow water on a plain sloping beach.
% Use HDG by hybridizing the conservative trace state.
pde.model = "ModelD";
pde.modelfile = "pdemodel";

% Compute on CPU with a single rank by default.
pde.platform = "cpu";
pde.mpiprocs = 1;
pde.hybrid = 1;
pde.Cxxpreprocessing = 0;

% Physical parameters from swe_1d.m / swe_1d_writeup.pdf.
g = 9.81;
L = 10000.0;
zbeach = 3.0;
zdeep = -7.0;
A_tide = 0.4;
T_tide = 1800.0;
omega = 2.0*pi/T_tide;
H_floor = 1e-3;
H_reg = 1e-3;
pde.physicsparam = [g, L, zbeach, zdeep, A_tide, omega, H_floor];

% Conservative settings are intentional because wet/dry fronts are stiff.
N = 400;
dt = 2.0;
nsteps = round(3.0*T_tide/dt);
%nsteps = 4000;

pde.porder = 2;
pde.torder = 2;
pde.nstage = 2;
pde.tau = 4;
pde.dt = dt*ones(1, nsteps);
pde.saveSolFreq = 10;
pde.soltime = pde.saveSolFreq:pde.saveSolFreq:nsteps;
pde.visdt = pde.saveSolFreq*dt;

pde.GMRESrestart = 30;
pde.linearsolvertol = 1e-8;
pde.linearsolveriter = 30;
pde.NLtol = 1e-8;
pde.NLiter = 3;

pde.AV = 1;
pde.frozenAVflag = 1;
pde.AVsmoothingIter = 2;

% Uniform 1D mesh on [0, L].
[mesh.p, mesh.t] = linemesh(N);
mesh.p = L*mesh.p;
mesh.boundaryexpr = {@(p) abs(p(1,:)) < 1e-8, @(p) abs(p(1,:) - L) < 1e-8};
mesh.boundarycondition = [1; 2];

mesh.f = facenumbering(mesh.p,mesh.t,pde.elemtype,mesh.boundaryexpr,mesh.periodicexpr);
dgnodes = createdgnodes(mesh.p, mesh.t, mesh.f, mesh.curvedboundary, ...
    mesh.curvedboundaryexpr, pde.porder);
x = dgnodes(:,1,:);
[H0, Hu0] = initsol(x, pde.physicsparam);
% mesh.udg = H0;
% mesh.udg(:,2,:) = Hu0;
mesh.vdg = 0*x;

% Generate and run the Exasim model.
[sol, pde, mesh, master, dmd] = exasim(pde, mesh);

% Simple postprocessing of the final state.
H = squeeze(sol(:,1,:,end));
q = squeeze(sol(:,2,:,end));

zb = localbathymetry(x, pde.physicsparam);
zb = squeeze(zb);

for i = 1:size(sol,4)
  H = squeeze(sol(:,1,:,i));
  q = squeeze(sol(:,2,:,i));
  Hx = squeeze(sol(:,3,:,i));
  qx = squeeze(sol(:,4,:,i));

  gamma = 1e3;
  Hp = H.*(atan(gamma*(H))/pi + 0.5) - atan(gamma)/pi + 0.5;    
  dHpdH = (atan(gamma*H)/pi + 0.5) + H.*(gamma./(pi*(1 + gamma^2*H.^2)));
  Hpx = dHpdH .* Hx;
  den = (Hp + H_floor);
  u = q./den;
  ux = (qx.*den - q.*Hpx) ./ (den.^2);

  eta = H + zb;  
  
  s = ux;
  s = s/0.001;
  s = limiting(s,0.0,0.9,gamma,0.0);
  s = s/0.9;
  S0 = 0.2; 
  av = (s-S0).*(atan(gamma*(s-S0))/pi + 0.5) - atan(gamma)/pi + 0.5;    

  figure(1); clf; plot(x(:), H(:));
  figure(2); clf; plot(x(:), q(:));
  figure(3); clf; plot(x(:), eta(:));
  figure(4); clf; plot(x(:), u(:));
  figure(5); clf; plot(x(:), av(:));
  pause(0.25);
end

% xplot = x(:);
% Hplot = H(:);
% etaplot = eta(:);
% uplot = u(:);
% zbplot = zb(:);
% [xplot, ind] = sort(xplot);
% Hplot = Hplot(ind);
% etaplot = etaplot(ind);
% uplot = uplot(ind);
% zbplot = zbplot(ind);
% 
% figure(1); clf;
% subplot(3,1,1);
% plot(xplot, Hplot, "LineWidth", 1.5);
% xlabel("x [m]");
% ylabel("H [m]");
% title("Water depth");
% grid on;
% 
% subplot(3,1,2);
% plot(xplot, etaplot, "b-", "LineWidth", 1.5); hold on;
% plot(xplot, zbplot, "k--", "LineWidth", 1.0);
% xlabel("x [m]");
% ylabel("elevation [m]");
% legend("free surface", "bed", "Location", "best");
% title("Free surface and bathymetry");
% grid on;
% 
% subplot(3,1,3);
% plot(xplot, uplot, "LineWidth", 1.5);
% xlabel("x [m]");
% ylabel("u [m/s]");
% title("Depth-averaged velocity");
% grid on;

disp("Done!");

function zb = localbathymetry(x, mu)
L = mu(2);
zbeach = mu(3);
zdeep = mu(4);

slope = (zdeep - zbeach)/L;
zb = zbeach + slope*x;
end
