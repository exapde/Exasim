% Specify an Exasim version to run
version = "Version0.2";

% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde1 structure and mesh1 structure
[pde1,mesh1] = initializeexasim(version);

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde1.model = "ModelD";          % ModelC, ModelD, ModelW
pde1.modelfile = "pdemodel1";    % name of a file defining the PDE model
pde1.modelnumber = 1;

% Choose computing platform and set number of processors
%pde1.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde1.mpiprocs = 1;              % number of MPI processors

% Set discretization parameters, physical parameters, and solver parameters
pde1.porder = 3;          % polynomial degree
pde1.physicsparam = 1;    % unit thermal conductivity
pde1.tau = 1.0;           % DG stabilization parameter

% create a grid of 8 by 8 on the unit square
[mesh1.p,mesh1.t] = squaremesh(8,8,1,1);
% expressions for domain boundaries
mesh1.boundaryexpr = {@(p) abs(p(2,:))<1e-8, @(p) abs(p(1,:)-1)<1e-8, @(p) abs(p(2,:)-1)<1e-8, @(p) abs(p(1,:))<1e-8};
mesh1.boundarycondition = [1;1;1;1]; % Set boundary condition for each boundary

% generate code for PDE model 1
[pde1,mesh1,master1,dmd1] = generatecode(pde1,mesh1);

