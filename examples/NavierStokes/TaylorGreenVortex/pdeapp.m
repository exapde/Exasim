% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
%pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 8;              % number of MPI processors
pde.hybrid = 1;

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 2;          % polynomial degree
pde.torder = 3;          % time-stepping order of accuracy
pde.nstage = 3;          % time-stepping number of stages
pde.dt = 0.16*ones(1,200);   % time step sizes
pde.soltime = 1:length(pde.dt); % steps at which solution are collected
pde.visdt = 0.05; % visualization timestep size

gam = 1.4;                      % specific heat ratio
Re = 1600;                      % Reynolds number
Pr = 0.71;                      % Prandtl number    
Minf = 0.2;                     % Mach number
pde.physicsparam = [gam Re Pr Minf];
pde.tau = 2.0;                  % DG stabilization parameter
pde.GMRESrestart=50;
pde.linearsolvertol=1e-6;
pde.linearsolveriter=51;
pde.preconditioner=1;
pde.precMatrixType=2;
pde.NLiter=2;
pde.ppdegree = 0;

% create a grid of 10 by 10 on the unit square
[mesh.p,mesh.t] = cubemesh(24,24,24,1);
mesh.p = 2*pi*mesh.p;
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(2,:))<1e-8, @(p) abs(p(1,:)-2*pi)<1e-8, @(p) abs(p(2,:)-2*pi)<1e-8, @(p) abs(p(1,:))<1e-8, @(p) abs(p(3,:))<1e-8, @(p) abs(p(3,:)-2*pi)<1e-8};
mesh.boundarycondition = [1;1;1;1;1;1];
% Set periodic boundary conditions
mesh.periodicexpr = {2, @(p) p([2 3],:), 4, @(p) p([2 3],:); 1, @(p) p([1 3],:), 3, @(p) p([1 3],:); 5, @(p) p([1 2],:), 6, @(p) p([1 2],:)};

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh] = exasim(pde,mesh);

% visualize the numerical solution of the PDE model using Paraview
% pde.visscalars = {"density", 1, "energy", 5};  % list of scalar fields for visualization
% pde.visvectors = {"momentum", [2, 3, 4]}; % list of vector fields for visualization
% xdg = vis(sol,pde,mesh); % visualize the numerical solution
% disp("Done!");


