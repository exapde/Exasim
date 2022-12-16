% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

pde.cpucompiler="/opt/homebrew/Cellar/llvm@12/12.0.1_1/bin/clang++";

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelC";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel_euler";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
%pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors

%% Set discretization parameters
pde.porder = 1;          % polynomial degree
pde.torder = 1;          % time-stepping order of accuracy
pde.nstage = 1;           % time-stepping number of stages
nt = 10000;
pde.dt = 5e-4*ones(1,nt);   % time step sizes
pde.saveSolFreq = nt/100;
pde.soltime = pde.saveSolFreq:pde.saveSolFreq:length(pde.dt); % steps at which solution are collected
pde.tau = 4.0;

%% Mesh
nDiv = 1024;

[mesh.p,mesh.t] = linemesh(nDiv-1);
a = 0; b = 1;
mesh.p = a + (b-a)*mesh.p;
% mesh.p = loginc(mesh.p, 5);
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(1,:)-a)<1e-16, @(p) abs(p(1,:) - b)<1e-8}; %TODO: double check boundaries
mesh.boundarycondition = [1;2]; % Set boundary condition for each boundary


%% Solver params
pde.linearsolveriter = 40;
pde.GMRESrestart = 20;
pde.matvectol = 1e-3;
pde.NLiter = 3;
pde.NLtol = 1e-11;
%%  Problem params

% Postshock flow conditions 
% Dimensional quantities used in nonequilibrium case
rho_post = 0.0085;
u_post = 542.0;
rhou_post = rho_post*u_post;
rhoE_post = 52157.2;

p_post = 16101;
p_outflow = 16894;
% pde.physicsparam = [rho_post(:)'/rho_inf,...
%                     rhou_post/(rho_inf*u_inf),...
%                     rhoE_post/(rho_inf*u_inf^2)];
% pout = 2.104849700000000e+04/(rho_inf*u_inf^2);

% Assumed equilibrium conditions used for outflow and initial conditions.  
pde.physicsparam = [1, 0, 16.8925]; 

rhoin = 1.0;
uin = 1.0;
rhoEin = 20.888;
pin = 0.4 * (rhoEin - 1/2 * rhoin * uin^2);

% pout = p_outflow/(rho_post*u_post*u_post);
pout = 2.0;
% pout = 3.0;
Minf = 0.444658;

Tin = 1.4 * Minf^2 * pin/rhoin;

pde.externalparam = [rho_post, 1.0, pout, Minf, Tin];

%% call exasim to generate and run C++ code to solve the PDE model
% [sol,pde,mesh,master,dmd] = exasim(pde,mesh);

% search compilers and set options
pde = setcompilers(pde);       

% generate input files and store them in datain folder
[pde,mesh,master,dmd] = preprocessing(pde,mesh);

% generate source codes and store them in app folder
gencode(pde);

compilerstr = compilecode(pde);

% run executable file to compute solution and store it in dataout folder
runstr = runcode(pde);

