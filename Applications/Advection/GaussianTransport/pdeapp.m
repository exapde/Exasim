% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelC";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
%pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 4;          % polynomial degree
pde.torder = 3;          % time-stepping order of accuracy
pde.nstage = 3;          % time-stepping number of stages
pde.physicsparam = [1 1];    % convective velocity
pde.tau = 1.0;               % DG stabilization parameter
pde.dt = 0.01*ones(1,200);   % time step sizes
pde.soltime = 1:length(pde.dt); % steps at which solution are collected
pde.visdt = 0.01; % visualization timestep size

% create a grid of 10 by 10 on the unit square
[mesh.p,mesh.t] = squaremesh(10,10,1,1);
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(2,:))<1e-8, @(p) abs(p(1,:)-1)<1e-8, @(p) abs(p(2,:)-1)<1e-8, @(p) abs(p(1,:))<1e-8};
mesh.boundarycondition = [1;1;1;1]; 
% Set periodic boundary conditions
mesh.periodicexpr = {2, @(p) p(2,:), 4, @(p) p(2,:); 1, @(p) p(1,:), 3, @(p) p(1,:)};

% % call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh] = exasim(pde,mesh);

% visualize the numerical solution of the PDE model using Paraview
pde.visscalars = {"temperature", 1};  % list of scalar fields for visualization
pde.visvectors = {}; % list of vector fields for visualization
xdg = vis(sol,pde,mesh); % visualize the numerical solution
disp("Done!");


