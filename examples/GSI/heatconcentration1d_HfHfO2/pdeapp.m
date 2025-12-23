% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";         % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";   % name of a file defining the PDE model

% Choose computing platform and set number of processors
pde.platform = "cpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;             % number of MPI processors
pde.hybrid = 1;               % 0 -> LDG, 1 -> HDG
pde.Cxxpreprocessing = 0;

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;               % polynomial degree
pde.pgauss = 2*pde.porder;    % Gauss quadrature degree
pde.torder = 1;               % time-stepping order of accuracy
pde.nstage = 1;               % time-stepping number of stages
pde.tau = 1.0;                % DG stabilization parameter

L = 5e-2;                % m physical length of the solid 
Twall = 300;             % K wall temperature  
T0  = 1000;              % K  temperature used for nondimensionalization
R = 8.314;               % J/(K * mol) ideal gas constant
kB = 8.617e-5;           % eV/K Boltzmann Constant
k = 1e-6;                % 1/s reaction rate
H0 = 1.415*1e10;         % J/m^3  chemistry modulus
mu0 = -104.9 * 1e3;      % J/mol reference chemical potential 
C0 = 7.45 * 1e4;         % mol/m^3 reference concentration
D0 = 7.8*1e-5;          % m^2/s diffusion rate constant
H = 1.46;                % eV diffusion barrier energy

rhoHf = 12500;           % kg/m^3 Hf density
cpHf  = 220;             % J/(kg * K) Hf heat capacity
kappaHf = 23;            % J/(s * m * K) Hf thermal conductivity
CvHf = cpHf*rhoHf;       % J/(m^3 * K) Hf volumetric specific heat;
alphaHf = kappaHf/CvHf;  % m^2/s Hf thermal diffusivity

rhoHfO2 = 9500;                % kg/m^3 HfO2 density
cpHfO2  = 120;                 % J/(kg * K) HfO2 heat capacity
kappaHfO2 = 2.0;               % J/(s * m * K) HfO2 thermal conductivity
CvHfO2 = cpHfO2*rhoHfO2;       % J/(m^3 * K) HfO2 volumetric specific heat;
alphaHfO2 = kappaHfO2/CvHfO2;  % m^2/s HfO2 thermal diffusivity

tc = L*L/alphaHf;              % s physical final time  
Dc = tc*D0/(L*L);              % dimensionless diffusion constant
DTHf = tc*alphaHf/(L*L);       % dimensionless thermal diffusivity for Hf
DTHfO2 = tc*alphaHfO2/(L*L);   % dimensionless thermal diffusivity for HfO2

heatflux = 1.32;         % dimensionless heat flux

% kappa * dT/dx = J/(s * m * K) * (K / m) = J/(m^2 * s) 
% kappa/(rho * cp) dT/dx = alpha * dT/dx   
% kappa dT/dx = rho * cp * (alpha * dT/dx)
heatflux_phy = heatflux * alphaHf * T0/L;         % K*m/s physical heat flux based on thermal diffusivity
heatfluxcond_phy = heatflux_phy * rhoHf * cpHf;   % J/(m^2 * s)  physical heat flux based on thermal conductivity

betaT = 1;
mu = [L Twall T0 R kB k H0 mu0 C0 D0 H tc Dc DTHf CvHf DTHfO2 CvHfO2 heatflux betaT];

pde.physicsparam = mu;   % phyiscal parameters
pde.ppdegree = 10;       % degree of polynomial preconditioner
pde.RBdim = 0;           % reduced basis dimension
pde.tau = 1;             % DG stabilization parameter
pde.NLiter = 3;

nsteps = 800;                   % number of time steps
pde.dt = 0.0001*(1:1:nsteps).^2; % timestep sizes
pde.soltime = 1:length(pde.dt); % steps at which solution are collected
pde.linearsolvertol = 1e-8;     % GMRES tolerance

% create a 1D grid from 0 to 1
[mesh.p,mesh.t] = linemesh(300);

% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(1))<1e-8, @(p) abs(p(1)-1<1e-8)};

mesh.boundarycondition = [1;2]; %set boundary condition for each boundary

pde.gencode = 1;
% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh,master,dmd] = exasim(pde,mesh);

figure(1); clf; hold on;
for i = 1:10:nsteps
  u = sol(:,1,:,i);
  plot(mesh.dgnodes(:),u(:));
end

figure(2); clf; hold on;
for i = 10:10:nsteps
  u = sol(:,2,:,i);
  plot(mesh.dgnodes(:),u(:));
end

return;

L0 = 5;
nsteps = length(pde.dt);

I = [10 20 30 100 160 250 400];
t=cumsum(pde.dt)*tc;
ti=t(I)/3600;

figure(1); clf; hold on;
for i = 1:length(I)
  u1 = sol(:,1,:,I(i));
  plot(L0*mesh.dgnodes(:),T0*u1(:),'-','LineWidth',2); 
end
set(gca,'FontSize',22); 
set(gca,'LooseInset',get(gca,'TightInset'))
xlabel("$\mbox{Solid domain (cm)}$", 'interpreter', 'latex', 'FontSize', 26);
ylabel("$T (\mbox{K})$", 'interpreter', 'latex', 'FontSize', 26);
set(gca, 'XTick', [0:1:10]/2);
set(gca, 'YTick', [0:1:10]*300);
set(gca,'TickLength',[0.02, 0.02])
axis([0 L0 0 3*T0])
box on;
grid on;
%leg = legend({'$t = 0.0032 \mbox{ } \mathrm{h}$','$t = 0.0238 \mbox{ } \mathrm{h}$','$t = 0.0785 \mbox{ } \mathrm{h}$','$t = 2.8094 \mbox{ } \mathrm{h}$','$t = 43.51 \mbox{ } \mathrm{h}$','$t = 177.8 \mbox{ } \mathrm{h}$','$t = 1420 \mbox{ } \mathrm{h}$'}, 'interpreter', 'latex', 'FontSize', 20, 'Location', 'NE');
leg = legend({'$t = 0.0032 \mbox{ } \mathrm{h}$','$t = 0.0238 \mbox{ } \mathrm{h}$','$t = 0.0785 \mbox{ } \mathrm{h}$','$t = 2.8094 \mbox{ } \mathrm{h}$','$t = 11.44 \mbox{ } \mathrm{h}$','$t = 43.51 \mbox{ } \mathrm{h}$','$t = 178.8 \mbox{ } \mathrm{h}$'}, 'interpreter', 'latex', 'FontSize', 20, 'Location', 'NE');
leg.ItemTokenSize = [40,10];
print -dpng oxidation1d_HfHfO2_temp_expdiff.png

I = [10 60 110 160 223 320 500];
ti=t(I)/3600;

figure(2); clf; hold on;
for i = 1:length(I)
  u1 = sol(:,2,:,I(i));
  plot(L0*mesh.dgnodes(:),u1(:),'-','LineWidth',2); 
end
set(gca,'FontSize',22); 
set(gca,'LooseInset',get(gca,'TightInset'))
xlabel("$\mbox{Solid domain (cm)}$", 'interpreter', 'latex', 'FontSize', 26);
ylabel("$\mathcal{E} = \mathcal{C}/\mathcal{C}_0$", 'interpreter', 'latex', 'FontSize', 26);
set(gca, 'XTick', [0:1:10]/2);
set(gca, 'YTick', [0:1:10]*0.1);
set(gca,'TickLength',[0.02, 0.02])
axis([0 L0 0 1])
box on;
grid on;
%leg = legend({'$t = 0.0032 \mbox{ } \mathrm{h}$','$t =2.8094 \mbox{ } \mathrm{h}$','$t = 19.134 \mbox{ } \mathrm{h}$','$t = 54.78 \mbox{ } \mathrm{h}$','$t = 129.6 \mbox{ } \mathrm{h}$','$t = 307.0 \mbox{ } \mathrm{h}$','$t = 1420 \mbox{ } \mathrm{h}$'}, 'interpreter', 'latex', 'FontSize', 20, 'Location', 'NE');
leg = legend({'$t = 0.0032 \mbox{ } \mathrm{h}$','$t = 0.6129 \mbox{ } \mathrm{h}$','$t = 3.7342 \mbox{ } \mathrm{h}$','$t = 11.443 \mbox{ } \mathrm{h}$','$t = 30.9 \mbox{ } \mathrm{h}$','$t = 91.12 \mbox{ } \mathrm{h}$','$t = 347 \mbox{ } \mathrm{h}$'}, 'interpreter', 'latex', 'FontSize', 20, 'Location', 'NE');
leg.ItemTokenSize = [40,10];
print -dpng oxidation1d_HfHfO2_volfrac_expdiff.png


