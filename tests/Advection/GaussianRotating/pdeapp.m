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
pde.dt = 0.025*ones(1,200);   % time step sizes
pde.soltime = 10:10:200; % steps at which solution are collected
pde.visdt = 0.025; % visualization timestep size

% create a grid of 10 by 10 on the unit square
[mesh.p,mesh.t] = squaremesh(20,20,1,1);
mesh.p = mesh.p - 0.5;
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(2,:)+0.5)<1e-4, @(p) abs(p(1,:)-0.5)<1e-4, @(p) abs(p(2,:)-0.5)<1e-4, @(p) abs(p(1,:)+0.5)<1e-4};
mesh.boundarycondition = [1;1;1;1]; 

% % call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh] = exasim(pde,mesh);

