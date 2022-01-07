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
pde.dt = 0.05*ones(1,200);   % time step sizes
pde.soltime = 10:10:200; % steps at which solution are collected
pde.visdt = 0.05; % visualization timestep size

pde.physicsparam = [5/3]; % specific heat ratio
pde.tau = 4;            % DG stabilization parameter

% create a grid of 16 by 16 on the unit square
[mesh.p,mesh.t] = squaremesh(16,16,1,1);
mesh.p = 10.0*mesh.p - 5.0;
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(2,:)+5)<1e-8, @(p) abs(p(1,:)-5)<1e-8, @(p) abs(p(2,:)-5)<1e-8, @(p) abs(p(1,:)+5)<1e-8};
mesh.boundarycondition = [1;1;1;1];
% Set periodic boundary conditions
mesh.periodicexpr = {2, @(p) p(2,:), 4, @(p) p(2,:); 1, @(p) p(1,:), 3, @(p) p(1,:)};

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh] = exasim(pde,mesh);

