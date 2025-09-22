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

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;               % polynomial degree
pde.pgauss = 2*pde.porder;    % Gauss quadrature degree
pde.torder = 1;               % time-stepping order of accuracy
pde.nstage = 1;               % time-stepping number of stages
pde.tau = 1.0;                % DG stabilization parameter

L = 5e-2;                % m physical length of the solid 
Twall = 500;             % K wall temperature  
T0  = 1000;              % K  temperature used for nondimensionalization
R = 8.314;               % J/(K * mol) ideal gas constant
kB = 8.617e-5;           % eV/K Boltzmann Constant
k = 1e-5;                % 1/s reaction rate
H0 = 1.415*1e10;         % J/m^3  chemistry modulus
mu0 = -104.9 * 1e3;      % J/mol reference chemical potential 
C0 = 7.45 * 1e4;         % mol/m^3 reference concentration
D0 = 9.28*1e-6;          % m^2/s diffusion rate constant
H = 1.46;                % eV diffusion barrier energy

rhoHf = 13000;           % kg/m^3 Hf density
cpHf  = 170;             % J/(kg * K) Hf heat capacity
kappaHf = 15;            % J/(s * m * K) Hf thermal conductivity
CvHf = cpHf*rhoHf;       % J/(m^3 * K) Hf volumetric specific heat;
alphaHf = kappaHf/CvHf;  % m^2/s Hf thermal diffusivity

tc = L*L/alphaHf;        % s physical final time  
Dc = tc*D0/(L*L);        % dimensionless diffusion constant
DTHf = tc*alphaHf/(L*L); % dimensionless thermal diffusivity

heatflux = 1.5;         % dimensionless heat flux

mu = [L Twall T0 R kB k H0 mu0 C0 D0 H tc Dc DTHf CvHf heatflux];

pde.physicsparam = mu;   % phyiscal parameters
pde.ppdegree = 10;       % degree of polynomial preconditioner
pde.RBdim = 0;           % reduced basis dimension
pde.tau = 1;             % DG stabilization parameter

nsteps = 200;                   % number of time steps
pde.dt = 0.001*(1:1:nsteps).^2; % timestep sizes
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
for i = 1:1:20
  u = sol(:,1,:,i);
  plot(mesh.dgnodes(:),u(:));
end

figure(2); clf; hold on;
for i = 10:10:nsteps
  u = sol(:,2,:,i);
  plot(mesh.dgnodes(:),u(:));
end


 

