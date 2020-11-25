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
pde.mpiprocs = 1;              % number of MPI processors

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;          % polynomial degree
pde.torder = 2;          % time-stepping order of accuracy
pde.nstage = 2;          % time-stepping number of stages
pde.dt = 0.01*ones(10000,1);   % time step sizes
pde.visdt = pde.dt(1);         % visualization timestep size
pde.saveSolFreq = 100;          % solution is saved every 100 time steps
pde.soltime = 20000; % steps at which solution are collected
pde.timestepOffset = 20000;
pde.paraview = "/Applications/ParaView-5.8.0.app/Contents/MacOS/paraview";

R0 = 10.0;
R1 = 11.5;

gam = 1.4;                      % specific heat ratio
Re = 500;                      % Reynolds number
Pr = 0.71;                      % Prandtl number    
Minf = 0.1;                     % Mach number
rbot = 1.0;                     % Density bottom boundary
Tbot = 1.0;                     % Temperature bottom surface
pbot = rbot*Tbot/(gam*Minf^2);  % Pressure bottom boundary
Ttop = 1.0;                     % Temperature top surface (equal to Tbot to start)
gravity = 0.7;                  % Normalized gravity acceleration
pde.physicsparam = [gam Re Pr Minf gravity rbot pbot Tbot Ttop R0 R1];
                   % 1  2  3   4      5      6    7    8    9  10 11
pde.tau = 2.0;                  % DG stabilization parameter
pde.GMRESrestart=30;            % number of GMRES restarts
pde.linearsolvertol=0.00001;     % GMRES tolerance
pde.linearsolveriter=31;        % number of GMRES iterations
pde.precMatrixType=2;           % preconditioning type!mv
pde.NLtol = 1e-8;               % Newton toleranccd dataoue
pde.NLiter = 3;                 % Newton iterations

% read a grid from a file
[mesh.p,mesh.t, mesh.dgnodes] = sphere_cube2(R0,R1,22,15,pde.porder);
% expressions for domain boundaries
%mesh.boundaryexpr = {@(p) abs(sqrt(sum(p.^2,1))-R0)<1e-8, @(p) abs(sqrt(sum(p.^2,1))-R1)<1e-8};
mesh.boundaryexpr = {@(p) abs(p(1,:).^2+p(2,:).^2+p(3,:).^2-R0^2)<1e-6, @(p) abs(p(1,:).^2+p(2,:).^2+p(3,:).^2-R1^2)<1e-6};
mesh.boundarycondition = [1 2];  % Inner, Outer
% expressions for curved boundaries
mesh.curvedboundary = [0 0];
mesh.curvedboundaryexpr = {@(p) sqrt(p(1,:).^2+p(2,:).^2+p(3,:).^2)-R0, @(p) sqrt(p(1,:).^2+p(2,:).^2+p(3,:).^2)-R1};

% % call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh] = exasim(pde,mesh);
% 
mesh.dgnodes = mesh.dgnodes + 10;

% visualize the numerical solution of the PDE model using Paraview
pde.visscalars = {"density", 1, "energy", 5};  % list of scalar fields for visualization
pde.visvectors = {"momentum", [2, 3, 4]}; % list of vector fields for visualization
xdg = vis(sol,pde,mesh); % visualize the numerical solution
disp("Done!");

