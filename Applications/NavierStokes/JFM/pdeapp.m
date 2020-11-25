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
pde.dt = 0.003*ones(1,10);   % time step sizes
pde.saveSolFreq = 10;          % solution is saved every 10 time steps
pde.soltime = 100:100:length(pde.dt); % steps at which solution are collected
pde.visdt = pde.dt(1);           % visualization timestep size

gam = 1.4;              % specific heat ratio
Minf = 0.3;             % Mach number
alpha = 3*pi/180;
rinf = 1.0;             % freestream density
uinf = cos(alpha);      % freestream x-horizontal velocity
vinf = 0.0;             % freestream y-horizontal velocity
winf = sin(alpha);             % freestream vertical velocity
pinf = 1/(gam*Minf^2);  % freestream pressure
rEinf = 0.5+pinf/(gam-1); % freestream energy
pde.physicsparam = [gam Minf rinf uinf vinf winf rEinf];
pde.tau = 2.0;          % DG stabilization parameter

pde.GMRESrestart=30;  % number of GMRES restarts
pde.linearsolvertol=0.0001; % GMRES tolerance
pde.linearsolveriter=31;  % number of GMRES iterations
pde.precMatrixType=2; % preconditioning type
pde.NLtol = 1e-7;  % Newton tolerance
pde.NLiter = 3;   % number of Newton iterations

% read a grid from a file
load("small_mesh/P");
load("small_mesh/T");
load("small_mesh/dgNodes");
mesh.p = P/1000;
mesh.t = T;
mesh.dgnodes = dgNodes/1000;
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) sqrt(sum(p.^2,1))<10, @(p) sqrt(sum(p.^2,1))>10};
mesh.boundarycondition = [1;2];

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh] = exasim(pde,mesh);

% visualize the numerical solution of the PDE model using Paraview
pde.visscalars = {"density", 1, "energy", 5};  % list of scalar fields for visualization
pde.visvectors = {"momentum", [2, 3, 4]}; % list of vector fields for visualization
xdg = vis(sol,pde,mesh); % visualize the numerical solution
disp("Done!");
