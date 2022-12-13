% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

pde.cpucompiler="/opt/homebrew/Cellar/llvm@12/12.0.1_1/bin/clang++";

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
pde.mpiprocs = 1;              % number of MPI processors

%% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;          % polynomial degree
pde.torder = 1;          % time-stepping order of accuracy
pde.nstage = 1;          % time-stepping number of stages

% tfinal = 0.000005; % 1 with 1e-4 doesn't work super well
tfinal = 0.0324;
% nt = 3240/10;
% dt = tfinal/nt;
% dt = 10*1e-5;
nt = 1000;
% nt = 150
% dt = tfinal / nt;
dt = 3.240000000000000e-06*10;
% pde.dt = [dt/100*ones(1,nt*10), dt*ones(1,nt*9/10)];
pde.dt = [dt*ones(1,nt)];
pde.saveSolFreq = ceil(nt/100);
% pde.saveSolFreq = 1;
pde.soltime = pde.saveSolFreq:pde.saveSolFreq:length(pde.dt); % steps at which solution are collected
% pde.timestepOffset = 10000;

% Solver params
pde.linearsolveriter = 40;
pde.GMRESrestart = 39;
pde.linearsolvertol = 1e-5;
pde.matvectol = 1e-7;
pde.NLiter = 5;
pde.RBdim = 5;
pde.NLtol = 1e-8;
pde.precMatrixType = 2;
pde.ptcMatrixType = 0;

%% Mutation information
% Mutation configuration 
pde.mutationflag = 1;
pde.mutationopts = {};
pde.mutationopts{1} = "air_5";
pde.mutationopts{2} = "ChemNonEq1T";
pde.mutationopts{3} = "RRHO";
pde.mutationpath = "/Users/rloekvh/Mutationpp";
nspecies = 5;

%% Stabilization and mesh
ndim = 1;
pde.tau = [5, 5, 5, 5, 5, 5, 5];
nDiv = 400;
[mesh.p,mesh.t] = linemesh(nDiv);
a = 0; b = 1;
mesh.p = a + (b-a)*mesh.p;

mesh.boundaryexpr = {@(p) abs(p(1,:)-a)<1e-16, @(p) abs(p(1,:) - b)<1e-8}; %TODO: double check boundaries
mesh.boundarycondition = [1;2]; % Set boundary condition for each boundary

%% Initial conditions and parameters
% TL=9000
Xi_inflow = [0,0,0,0,0];
u_inflow = 0.0; 
T_inflow = 9000.0;
p_outflow = 1e4;
pde.physicsparam = [
                        2.79122095e-02      % .
                        8.94162529e-03      % .
                        3.49305942e-05      % rho_i_L
                        1.58251700e-03      % . 
                        5.01765831e-07      % .
                        0.0                 % rhou_L
                        1425938.1573890178  % rhoE_L
                        8.83171348e-81      % .
                        4.62782812e-42      % .
                        3.11978853e-17      % rho_i_R
                        8.87231621e-02      % . 
                        2.69399686e-02      % .
                        0.0                 % rhou_R
                        -9783.724551855581  % rhoE_R
                        pde.porder          % p for h/p scaling
                        1.5               % bulk viscosity scaling parameter
                        0.0002             % cutoff dilitation
                        2.5                % maximum dilitation 18
                        Xi_inflow(:)     % 19 20 21 22 23
                        u_inflow         % 24
                        T_inflow         % 25
                        p_outflow        % 26
                    ];
% T = 4000;
% pde.physicsparam = [
%                         8.79041126e-05
%                         2.18935322e-02
%                         8.84111098e-03
%                         1.10501923e-01
%                         8.22504485e-03
%                         0.0                 % rhou_L
%                         863520.4456544879  % rhoE_L
%                         8.83171348e-81      % .
%                         4.62782812e-42      % .
%                         3.11978853e-17      % rho_i_R
%                         8.87231621e-02      % . 
%                         2.69399686e-02      % .
%                         0.0                 % rhou_R
%                         -9783.724551855581  % rhoE_R
%                         pde.porder          % p for h/p scaling
%                         1.5                % bulk viscosity scaling parameter
%                         0.02             % cutoff dilitation
%                         2.5                % maximum dilitation
%                         Xi_inflow(:)     % 19 20 21 22 23
%                         u_inflow         % 24
%                         T_inflow         % 25
%                         p_outflow        % 26
%                     ];


pde.externalparam = zeros(3*nspecies + 3 + ndim);
pde.externalparam(1) = 0.11566313067881752; % rho_inf
pde.externalparam(2) = 347.764810747843;    % u_inf
pde.externalparam(3) = 13988.341078772399;  % rhoE_inf
pde.externalparam(4) = 300;
pde.externalparam(5) = 1.9423965904207942e-05;
pde.externalparam(6) = 0.028469303787820466;

% mesh size
mesh.dgnodes = createnodes(mesh.p, mesh.t, pde.porder);
mesh.vdg = zeros(size(mesh.dgnodes,1),2,size(mesh.dgnodes,3));
mesh.vdg(:,2,:) = (1.0 / (nDiv)) * ones(size(mesh.dgnodes,1), 1, size(mesh.dgnodes,3));

pde.AV = 1;      
pde.nco = 10;
%% Call Exasim
% search compilers and set options
pde = setcompilers(pde);       

% generate input files and store them in datain folder
[pde,mesh,master,dmd] = preprocessing(pde,mesh);

% generate source codes and store them in app folder
gencode(pde);
cd("app");
eval('!./fixargsMPP_Mac.sh');
eval('!cp opuApp_MPP.cpp opuApp.cpp');
eval('!cp opuFlux_MPP.cpp opuFlux.cpp'); %X manually switch these 
eval('!cp opuFbou_MPP.cpp opuFbou.cpp'); %X
eval('!cp opuSource_MPP.cpp opuSource.cpp'); %X
eval('!cp opuOutput_MPP.cpp opuOutput.cpp'); %X
% eval('!cp opuOutput_viscTerms.cpp opuOutput.cpp');
eval('!cp opuAvfield_MPP.cpp opuAvfield.cpp') %X
eval('!cp opuInitu_MPP.cpp opuInitu.cpp') %X
% eval('!cp opuUbou_MPP.cpp opuUbou.cpp')
cd("..");

% compile source codes to build an executable file and store it in app folder
compilerstr = compilecode(pde);

% run executable file to compute solution and store it in dataout folder
runstr = runcode(pde);   
