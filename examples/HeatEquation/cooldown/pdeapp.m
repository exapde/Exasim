% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";       % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel"; % name of a file defining the PDE model

% Choose computing platform and set number of processors
pde.platform = "cpu";      % choose this option if you want to run the C++ code on Nvidia GPUs
pde.mpiprocs = 4;           % number of MPI processors
pde.hybrid = 1;
pde.debugmode = 0;

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;          % polynomial degree
pde.pgauss = 2*pde.porder;  % gauss quad order
pde.torder = 1;          % time-stepping order of accuracy
pde.nstage = 1;          % time-stepping number of stages
pde.physicsparam = 1;    % unit thermal conductivity
pde.tau = 1.0;           % DG stabilization parameter
pde.dt = 0.025*ones(1,20);   % time step sizes
pde.soltime = 1:length(pde.dt); % steps at which solution are collected
pde.visdt = 0.025; % visualization timestep size
pde.RBdim = 2;
pde.linearsolvertol = 1e-6; % GMRES tolerance

% create a grid of 8 by 8 on the unit square
[mesh.p,mesh.t] = squaremesh(8,8,1,1);
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(2,:))<1e-8, @(p) abs(p(1,:)-1)<1e-8, @(p) abs(p(2,:)-1)<1e-8, @(p) abs(p(1,:))<1e-8};
mesh.boundarycondition = [1;1;1;1]; % Set boundary condition for each boundary

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh,master,dmd] = exasim(pde,mesh);

% % visualize the numerical solution of the PDE model using Paraview
% pde.visscalars = {"temperature", 1};  % list of scalar fields for visualization
% pde.visvectors = {"temperature gradient", [2 3]}; % list of vector fields for visualization
% %xdg = vis(sol,pde,mesh); % visualize the numerical solution
% disp("Done!");

        
