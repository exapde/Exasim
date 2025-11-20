% initialize pde{1} structure and mesh{1} structure
[pde{1},mesh{1}] = initializeexasim();
pde{1}.buildpath = string(pwd);

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde{1}.model = "ModelD";          % ModelC, ModelD, ModelW
pde{1}.modelfile = "pdemodel1";    % name of a file defining the PDE model
%pde{1}.modelnumber = 1;
pde{1}.Cxxpreprocessing = 0;

% Choose computing platform and set number of processors
%pde{1}.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde{1}.mpiprocs = 1;              % number of MPI processors
pde{1}.hybrid = 1;
pde{1}.linearsolvertol = 1e-6;

% Set discretization parameters, physical parameters, and solver parameters
pde{1}.porder = 3;          % polynomial degree
pde{1}.physicsparam = 1;    % unit thermal conductivity
pde{1}.tau = 1.0;           % DG stabilization parameter

% create a grid of 8 by 8 on the unit square
[mesh{1}.p,mesh{1}.t] = squaremesh(4,8,1,1);
mesh{1}.p(1,:) = 0.5*mesh{1}.p(1,:);

% expressions for domain boundaries
mesh{1}.boundaryexpr = {@(p) abs(p(2,:))<1e-8, @(p) abs(p(1,:)-0.5)<1e-8, @(p) abs(p(2,:)-1)<1e-8, @(p) abs(p(1,:))<1e-8};
mesh{1}.boundarycondition = [1;2;1;1]; % Set boundary condition for each boundary

% % call exasim to generate and run C++ code to solve the PDE models
% [sol1,pde1,mesh1,master1,dmd1,compilerstr1,runstr1] = exasim(pde{1},mesh{1});
% 
% mesh1.porder = pde1.porder;
% figure(1); clf; scaplot(mesh1,sol1(:,1,:),[],2,1); axis on; axis equal; axis tight;
% 

