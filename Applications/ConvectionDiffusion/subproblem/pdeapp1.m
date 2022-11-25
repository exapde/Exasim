
% initialize pde{1} structure and mesh{1} structure
[pde{1},mesh{1}] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde{1}.model = "ModelD";          % ModelC, ModelD, ModelW
pde{1}.modelfile = "pdemodel1";    % name of a file defining the PDE model
pde{1}.modelnumber = 1;

% Choose computing platform and set number of processors
%pde{1}.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde{1}.mpiprocs = 1;              % number of MPI processors

% Set discretization parameters, physical parameters, and solver parameters
pde{1}.porder = 3;          % polynomial degree
pde{1}.physicsparam = 0.1;    % unit thermal conductivity
pde{1}.tau = 1.0;           % DG stabilization parameter

pde{1}.NLtol = 1.0000e-06;
pde{1}.linearsolvertol = 1.0000e-04;       
pde{1}.ppdegree = 20;
pde{1}.precMatrixType = 2;

% set indices to obtain v from the solutions of the other PDE models 
% first column : model index
% second column: solution index
pde{1}.vindx = [2 1; 2 2; 2 3];
pde{1}.subproblem = 1;

% create a grid of 8 by 8 on the unit square
[mesh{1}.p,mesh{1}.t] = squaremesh(32,32,1,1);
% expressions for domain boundaries
mesh{1}.boundaryexpr = {@(p) abs(p(2,:))<1e-8, @(p) abs(p(1,:)-1)<1e-8, @(p) abs(p(2,:)-1)<1e-8, @(p) abs(p(1,:))<1e-8};
mesh{1}.boundarycondition = [1;1;1;1]; % Set boundary condition for each boundary



