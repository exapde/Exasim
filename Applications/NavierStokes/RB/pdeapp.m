% Specify an Exasim version to run
version = "Version0.1";

% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim(version);

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
%pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 4;              % number of MPI processors

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;          % polynomial degree
pde.torder = 2;          % time-stepping order of accuracy
pde.nstage = 2;          % time-stepping number of stages
pde.dt = 0.001*ones(200,1);   % time step sizes
pde.visdt = pde.dt(1);         % visualization timestep size
pde.saveSolFreq = 10:10:200;          % solution is saved every 100 time steps
pde.soltime = 10:10:200; % steps at which solution are collected
pde.timestepOffset = 0;

gam = 1.4;                      % specific heat ratio
Re = 600;                      % Reynolds number
Pr = 0.71;                      % Prandtl number    
Minf = 0.1;                     % Mach number
rbot = 1.0;                     % Density bottom boundary
Tbot = 1.0;                     % Temperature bottom surface
pbot = rbot*Tbot/(gam*Minf^2);  % Pressure bottom boundary
Ttop = 0.5;                     % Temperature top surface (equal to Tbot to start)
gravity = 0.7;                  % Normalized gravity acceleration
pde.physicsparam = [gam Re Pr Minf gravity rbot pbot Tbot Ttop];
pde.tau = 5.0;                  % DG stabilization parameter
pde.GMRESrestart=30;            % number of GMRES restarts
pde.linearsolvertol=0.0001;     % GMRES tolerance
pde.linearsolveriter=31;        % number of GMRES iterations
pde.precMatrixType=2;           % preconditioning type
pde.NLtol = 1e-7;               % Newton toleranccd dataoue
pde.NLiter = 3;                 % Newton iterations

% read a grid from a file
[mesh.p,mesh.t] = squaremesh(15,15,0,1);
L = 5.0;
mesh.p(1,:) = L*mesh.p(1,:);
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(1,:))<1e-8, @(p) abs(p(2,:))<1e-8, @(p) abs(p(1,:)-L)<1e-8, @(p) abs(p(2,:)-1.0)<1e-8};
mesh.boundarycondition = [1;1;1;2];  % Left, Bottom, Right, Top
% Set periodic boundary conditions
mesh.periodicexpr = {1, @(p) p(2,:), 3, @(p) p(2,:)};

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh] = exasim(pde,mesh);

% visualize the numerical solution of the PDE model using Paraview
pde.visscalars = {"density", 1, "energy", 4};  % list of scalar fields for visualization
pde.visvectors = {"momentum", [2, 3]}; % list of vector fields for visualization
xdg = vis(sol,pde,mesh); % visualize the numerical solution
disp("Done!");

