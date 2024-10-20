% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Install/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelC";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
%pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors
pde.hybrid = 0;

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 4;          % polynomial degree
pde.torder = 3;          % time-stepping order of accuracy
pde.nstage = 3;          % time-stepping number of stages
pde.dt = 0.05*ones(1,200);   % time step sizes
pde.soltime = 1:length(pde.dt); % steps at which solution are collected
pde.visdt = 0.05; % visualization timestep size

gam = 1.4;                      % specific heat ratio
M_ref = sqrt(1/gam);            % Mach number
pde.physicsparam = [gam M_ref];
pde.tau = 1+1/M_ref;            % DG stabilization parameter

% create a grid of 10 by 10 on the unit square
[mesh.p,mesh.t] = squaremesh(10,10,1,1);
mesh.p = 10.0*mesh.p - 5.0;
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(2,:)+5)<1e-8, @(p) abs(p(1,:)-5)<1e-8, @(p) abs(p(2,:)-5)<1e-8, @(p) abs(p(1,:)+5)<1e-8};
mesh.boundarycondition = [1;1;1;1];
% Set periodic boundary conditions
mesh.periodicexpr = {2, @(p) p(2,:), 4, @(p) p(2,:); 1, @(p) p(1,:), 3, @(p) p(1,:)};

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh] = exasim(pde,mesh);

% visualize the numerical solution of the PDE model using Paraview
% pde.visscalars = {"density", 1, "energy", 4};  % list of scalar fields for visualization
% pde.visvectors = {"momentum", [2, 3]}; % list of vector fields for visualization
% xdg = vis(sol,pde,mesh); % visualize the numerical solution
% disp("Done!");
% 

mesh.porder = pde.porder;
mesh.dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);
for i = 1:size(sol,4)
  figure(1); clf; scaplot(mesh,sol(:,1,:,i),[-1 1],2,1); axis on; axis equal; axis tight;
end

