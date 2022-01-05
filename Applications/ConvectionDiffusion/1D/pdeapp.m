% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
%pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;          % polynomial degree
pde.physicsparam = [1 2*pi];    % unit thermal conductivity
pde.tau = 2*pi;           % DG stabilization parameter

% create a grid of 8 by 8 on the unit square
nDiv = 17;
%mesh.p = linspace(0,1,nDiv);
%mesh.t = [(1:nDiv-1); (2:nDiv)];
[mesh.p,mesh.t] = linemesh(nDiv-1);
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(1,:)-1)<1e-8, @(p) abs(p(1,:))<1e-8};
mesh.boundarycondition = [1;1]; % Set boundary condition for each boundary

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh,master,dmd] = exasim(pde,mesh);

% plot solution
u = sol(:,1,:);
dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);    
figure(1);clf;plot(dgnodes(:),u(:));


% % visualize the numerical solution of the PDE model using Paraview
% pde.visscalars = {"temperature", 1, "temperature gradient", 2};  % list of scalar fields for visualization
% pde.visvectors = {}; % list of vector fields for visualization
% xdg = vis(sol,pde,mesh); % visualize the numerical solution
% disp("Done!");

        
