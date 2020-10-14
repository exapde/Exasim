% Specify an Exasim version to run
version = "Version0.2";

% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde2 structure and mesh2 structure
[pde2,mesh2] = initializeexasim(version);

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde2.model = "ModelD";          % ModelC, ModelD, ModelW
pde2.modelfile = "pdemodel2";    % name of a file defining the PDE model
pde2.modelnumber = 2;

% Choose computing platform and set number of processors
%pde2.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde2.mpiprocs = 1;              % number of MPI processors

% Set discretization parameters, physical parameters, and solver parameters
pde2.porder = 3;          % polynomial degree
pde2.physicsparam = 1;    % unit thermal conductivity
pde2.tau = 1.0;           % DG stabilization parameter

% create a grid of 8 by 8 on the unit square
[mesh2.p,mesh2.t] = squaremesh(8,8,1,1);
% expressions for domain boundaries
mesh2.boundaryexpr = {@(p) abs(p(2,:))<1e-8, @(p) abs(p(1,:)-1)<1e-8, @(p) abs(p(2,:)-1)<1e-8, @(p) abs(p(1,:))<1e-8};
mesh2.boundarycondition = [1;1;1;1]; % Set boundary condition for each boundary

% generate code for PDE model 2
[pde2,mesh2,master2,dmd2] = generatecode(pde2,mesh2);

