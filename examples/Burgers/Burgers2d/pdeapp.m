% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Install/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
%pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors
pde.hybrid = 1;

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;          % polynomial degree
pde.physicsparam = [0.01 0 1 0.5 0];    % unit thermal conductivity
pde.tau = 1.0;           % DG stabilization parameter
pde.linearsolvertol = 1e-4;
pde.ppdegree = 10;

% create a grid of 8 by 8 on the unit square
[mesh.p,mesh.t] = squaremesh(32,32,1,1);
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(2,:))<1e-8, @(p) abs(p(1,:)-1)<1e-8, @(p) abs(p(2,:)-1)<1e-8, @(p) abs(p(1,:))<1e-8};
mesh.boundarycondition = [1;1;2;1]; % Set boundary condition for each boundary

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh] = exasim(pde,mesh);

% % visualize the numerical solution of the PDE model using Paraview
% pde.visscalars = {"temperature", 1};  % list of scalar fields for visualization
% pde.visvectors = {"temperature gradient", [2 3]}; % list of vector fields for visualization
% xdg = vis(sol,pde,mesh); % visualize the numerical solution
% disp("Done!");

mesh.porder = pde.porder;
mesh.dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);
figure(1); clf; scaplot(mesh,sol(:,1,:),[-1 1],2,1); axis on; axis equal; axis tight;

        
        
