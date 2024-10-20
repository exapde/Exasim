% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Install/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW 
pde.modelfile = "pdemodel4";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
%pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 4;          % polynomial degree
pde.torder = 2;          % time-stepping order of accuracy
pde.nstage = 2;          % time-stepping number of stages
pde.dt = 0.02*ones(1,5e3)/10;   % time step sizes
pde.visdt = 1; % visualization timestep size
pde.saveSolFreq = 50;          % solution is saved every 10 time steps
pde.soltime = 50:50:length(pde.dt); % steps at which solution are collected

gam = 1;   % gravity
pde.physicsparam = [gam];
pde.tau = 1;  % DG stabilization parameter
pde.GMRESrestart=15;            % number of GMRES restarts
pde.linearsolvertol=1e-12;       % GMRES tolerance
pde.linearsolveriter=16;        % number of GMRES iterations
pde.precMatrixType=2;           % preconditioning type
pde.NLtol = 1e-12;              % Newton tolerance
pde.NLiter=2;                   % Newton iterations

% parameters for differential algebraic equations
pde.dae_alpha = 1e-4;
pde.dae_beta = 0.0;

% create a grid of 10 by 10 on the unit square
[mesh.p,mesh.t] = squaremesh(128/2,128/2,1,1);
mesh.p = 4*pi*mesh.p - 2*pi;
% expressions for domain boundaries: bottom, right, top, left
mesh.boundaryexpr = {@(p) abs(p(2,:)+2*pi)<1e-8, @(p) abs(p(1,:)-2*pi)<1e-8, @(p) abs(p(2,:)-2*pi)<1e-8, @(p) abs(p(1,:)+2*pi)<1e-8};
%mesh.boundarycondition = [2;1;2;1]; % wall, perioidic, wall, periodic 
% Set periodic boundary conditions
%mesh.periodicexpr = {2, @(p) p(2,:), 4, @(p) p(2,:)};
mesh.boundarycondition = [1;1;1;1]; % wall, perioidic, wall, periodic 
mesh.periodicexpr = {2, @(p) p(2,:), 4, @(p) p(2,:); 1, @(p) p(1,:), 3, @(p) p(1,:)};

% search compilers and set options
pde = setcompilers(pde);       

% generate input files and store them in datain folder
[pde,mesh,master,dmd] = preprocessing(pde,mesh);

% generate source codes and store them in app folder
gencode(pde); 

% compile source codes to build an executable file and store it in app folder
compilerstr = compilecode(pde);

% % run executable file to compute solution and store it in dataout folder
runstr = runcode(pde);

sol = fetchsolution(pde,master,dmd);
sol(:,7,:,:) = sol(:,5,:,:)-sol(:,4,:,:);
pde.visscalars = {"vorticity", 7};
pde.visvectors = {"velocity", [1, 2]};
vis(sol,pde,mesh);

