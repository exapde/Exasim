% Specify an Exasim version to run
version = "Version0.1";

% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim(version);

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelW";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
%pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;          % polynomial degree
pde.physicsparam = 1;    % unit thermal conductivity
pde.tau = 1.0;           % DG stabilization parameter

% read a grid from a file
[mesh.p,mesh.t] = readmesh('grid.bin',0);
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(2,:)+12)<1e-3, @(p) abs(p(1,:)-12)<1e-3, @(p) abs(p(2,:)-12)<1e-3, @(p) abs(p(1,:)+12)<1e-3, @(p) abs(p(1,:))<20};
mesh.boundarycondition = [1;1;1;1;2]; % Set boundary condition for each boundary
% expressions for curved boundaries
mesh.curvedboundary = [0 0 0 0 1];
mesh.curvedboundaryexpr = {@(p) 0, @(p) 0, @(p) 0, @(p) 0, @(p) sqrt(p(1,:).*p(1,:)+p(2,:).*p(2,:))-1};
% call exasim to generate and run C++ code to solve the PDE model
% [sol,pde,mesh] = exasim(pde,mesh);

pde.elemtype = 0;
mesh.f = facenumbering(mesh.p,mesh.t,pde.elemtype,mesh.boundaryexpr,mesh.periodicexpr);
mesh.dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);    

meshplot(mesh.p,mesh.t);
axis equal; axis tight;
hold on;
x=mesh.dgnodes(:,1,:);
y=mesh.dgnodes(:,2,:);
plot(x(:),y(:),'*');

% % visualize the numerical solution of the PDE model using Paraview
% pde.visscalars = {"temperature", 1};  % list of scalar fields for visualization
% pde.visvectors = {"temperature gradient", [2 3]}; % list of vector fields for visualization
% vis(sol,pde,mesh); % visualize the numerical solution
% disp("Done!");
% 
