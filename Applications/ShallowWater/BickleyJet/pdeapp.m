% Specify an Exasim version to run
version = "Version0.1";

% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim(version);

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelC";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel2";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 4;          % polynomial degree
pde.torder = 2;          % time-stepping order of accuracy
pde.nstage = 2;          % time-stepping number of stages
pde.dt = 0.02*ones(1,5e3);   % time step sizes
pde.visdt = 1; % visualization timestep size
pde.saveSolFreq = 50;          % solution is saved every 10 time steps
pde.soltime = 50:50:length(pde.dt); % steps at which solution are collected

gam = 10;   % gravity
pde.physicsparam = [gam];
pde.tau = 1;  % DG stabilization parameter
pde.GMRESrestart=15;            % number of GMRES restarts
pde.linearsolvertol=1e-12;       % GMRES tolerance
pde.linearsolveriter=16;        % number of GMRES iterations
pde.precMatrixType=2;           % preconditioning type
pde.NLtol = 1e-12;              % Newton tolerance
pde.NLiter=2;                   % Newton iterations

% create a grid of 10 by 10 on the unit square
[mesh.p,mesh.t] = squaremesh(128/2,128/2,1,1);
mesh.p = 4*pi*mesh.p - 2*pi;
% expressions for domain boundaries: bottom, right, top, left
mesh.boundaryexpr = {@(p) abs(p(2,:)+2*pi)<1e-8, @(p) abs(p(1,:)-2*pi)<1e-8, @(p) abs(p(2,:)-2*pi)<1e-8, @(p) abs(p(1,:)+2*pi)<1e-8};
% Set periodic boundary conditions
mesh.boundarycondition = [1;1;1;1]; 
mesh.periodicexpr = {2, @(p) p(2,:), 4, @(p) p(2,:); 1, @(p) p(1,:), 3, @(p) p(1,:)};

% % call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh] = exasim(pde,mesh);
 
% % visualize the numerical solution of the PDE model using Paraview
pde.visscalars = {"density", 1};  % list of scalar fields for visualization
pde.visvectors = {"velocity", [2, 3]}; % list of vector fields for visualization
xdg = vis(sol,pde,mesh); % visualize the numerical solution
disp("Done!");


