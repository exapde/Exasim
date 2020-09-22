% Specify an Exasim version to run
version = "Version0.1";

% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim(version);

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelC";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
%pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 2;          % polynomial degree
pde.torder = 1;          % time-stepping order of accuracy
pde.nstage = 1;          % time-stepping number of stages
pde.dt = [0.01 0.1 1 10]/10;   % time step sizes
pde.soltime = 1:length(pde.dt); % steps at which solution are collected
pde.visdt = 1;           % visualization timestep size
pde.GMRESrestart = 100;
pde.precMatrixType = 2;

gam = 1.4;              % specific heat ratio
Minf = 0.4;             % Mach number
rinf = 1.0;             % freestream density
uinf = 1.0;             % freestream horizontal velocity
vinf = 0.0;             % freestream vertical velocity
pinf = 1/(gam*Minf^2);  % freestream pressure
rEinf = 0.5+pinf/(gam-1); % freestream energy
pde.physicsparam = [gam Minf rinf uinf vinf rEinf];
pde.tau = 2.0;          % DG stabilization parameter

% read a grid from a file
[mesh.p,mesh.t] = readmesh('grid.bin',0);
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) sqrt((p(1,:)-.5).^2+p(2,:).^2)<3, @(p) abs(p(1,:))<20};
mesh.boundarycondition = [1;2];
% expressions for curved boundaries
mesh.curvedboundary = [1 0];
mesh.curvedboundaryexpr = {@(p) p(2,:).^2-(5*0.01*12*(0.29690*sqrt(abs(p(1,:)))-0.12600*p(1,:)-0.35160*p(1,:).^2+0.28430*p(1,:).^3-0.10150*p(1,:).^4)).^2, @(p) 0};

% mesh.f = facenumbering(mesh.p,mesh.t,1,mesh.boundaryexpr,mesh.periodicexpr);
% [dgnodes,elemtype,perm] = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);
% meshplot(mesh.p,mesh.t); axis equal; axis tight;
% x = dgnodes(:,1,:); y = dgnodes(:,2,:);
% hold on; plot(x(:), y(:), 'o');
% x=linspace(0,1,10000);
% y=5*0.01*12*(0.29690*sqrt(x)-0.12600*x-0.35160*x.^2+0.28430*x.^3-0.10150*x.^4);
% plot(x(:), y(:), '-r');
% plot(x(:), -y(:), '-r');

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh] = exasim(pde,mesh);

% visualize the numerical solution of the PDE model using Paraview
pde.visscalars = {"density", 1, "energy", 4};  % list of scalar fields for visualization
pde.visvectors = {"momentum", [2, 3]}; % list of vector fields for visualization
xdg = vis(sol,pde,mesh); % visualize the numerical solution
disp("Done!");
