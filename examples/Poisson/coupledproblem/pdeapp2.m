% initialize pde{2} structure and mesh{2} structure
[pde{2},mesh{2}] = initializeexasim();
pde{2}.buildpath = string(pwd);

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde{2}.model = "ModelD";          % ModelC, ModelD, ModelW
pde{2}.modelfile = "pdemodel2";    % name of a file defining the PDE model
%pde{2}.modelnumber = 0;
pde{2}.Cxxpreprocessing = 0;

% Choose computing platform and set number of processors
%pde{2}.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde{2}.mpiprocs = 1;              % number of MPI processors
pde{2}.hybrid = 1;
pde{2}.linearsolvertol = 1e-6;

% Set discretization parameters, physical parameters, and solver parameters
pde{2}.porder = 3;          % polynomial degree
pde{2}.physicsparam = 1;    % unit thermal conductivity
pde{2}.tau = 1.0;           % DG stabilization parameter

% create a grid of 8 by 8 on the unit square
[mesh{2}.p,mesh{2}.t] = squaremesh(4,8,1,1);
mesh{2}.p(1,:) = 0.5*mesh{2}.p(1,:) + 0.5;

% expressions for domain boundaries
mesh{2}.boundaryexpr = {@(p) abs(p(2,:))<1e-8, @(p) abs(p(1,:)-1)<1e-8, @(p) abs(p(2,:)-1)<1e-8, @(p) abs(p(1,:)-0.5)<1e-8};
mesh{2}.boundarycondition = [1;1;1;2]; % Set boundary condition for each boundary

% % call exasim to generate and run C++ code to solve the PDE models
% [sol2,pde2,mesh2,master2] = exasim(pde{2},mesh{2});
% 
% mesh2.porder = pde2.porder;
% figure(2); clf; scaplot(mesh2,sol2(:,1,:),[],2,1); axis on; axis equal; axis tight;


