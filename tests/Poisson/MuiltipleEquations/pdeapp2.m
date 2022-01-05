% initialize pde{2} structure and mesh{2} structure
[pde{2},mesh{2}] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde{2}.model = "ModelD";          % ModelC, ModelD, ModelW
pde{2}.modelfile = "pdemodel2";    % name of a file defining the PDE model
pde{2}.modelnumber = 2;

% Choose computing platform and set number of processors
%pde{2}.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde{2}.mpiprocs = 1;              % number of MPI processors

% Set discretization parameters, physical parameters, and solver parameters
pde{2}.porder = 3;          % polynomial degree
pde{2}.physicsparam = 1;    % unit thermal conductivity
pde{2}.tau = 1.0;           % DG stabilization parameter

% set indices to obtain v from the solutions of the other PDE models 
% first column : model index
% second column: solution index
pde{2}.vindx = [1 1];
%pde{2}.vindx = [1 1; 1 2; 1 3];

% create a grid of 8 by 8 on the unit square
[mesh{2}.p,mesh{2}.t] = squaremesh(8,8,1,1);
% expressions for domain boundaries
mesh{2}.boundaryexpr = {@(p) abs(p(2,:))<1e-8, @(p) abs(p(1,:)-1)<1e-8, @(p) abs(p(2,:)-1)<1e-8, @(p) abs(p(1,:))<1e-8};
mesh{2}.boundarycondition = [1;1;1;1]; % Set boundary condition for each boundary


