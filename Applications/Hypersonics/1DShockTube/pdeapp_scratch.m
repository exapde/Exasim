% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

pde.cpucompiler="/opt/homebrew/Cellar/llvm@12/12.0.1_1/bin/clang++";

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelC";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel_inflow";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
pde.mpiprocs = 1;              % number of MPI processors

%% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 1;          % polynomial degree
pde.torder = 1;          % time-stepping order of accuracy
pde.nstage = 1;          % time-stepping number of stages
% time history
% 1000 to 0.01
% 1000 to 0.01
% 1000 to 0.01...no these don't work super well -__-
tfinal = 0.1; % 1 with 1e-4 doesn't work super well
% nt = 1000;
% dt = tfinal/nt;
% dt = 1e-4
% pde.dt = [dt/100*ones(1,nt/2), dt*ones(1,nt/2)];%, dt*2*ones(1,nt)];   % time step sizes
pde.dt = [(1e-3)*ones(1,20e4)];
% time history:
%  5e-6 x 5000 steps
% 5e-4 x 5000 steps
% Total so far: 2.525
% 5e-4 x 10000 steps
% 7.525
% 5e-3 x 10000 steps NOT GREAT: Oscillations and not stabilized yet
% 5e-2 x 4000 this does not seem right!
% 1e-2 x 1000
% 5e-2 x 1000
% 1e-3 x 10e3
% let's just keep going
% 1e-3 x ~85000
% pde.dt = [dt/100*ones(1,nt*10), dt*ones(1,nt*9/10)];
% pde.dt = [5*dt*ones(1,nt)]
% pde.saveSolFreq = nt/100;
pde.saveSolFreq = 500;
% pde.saveSolFreq = 1;
% pde.saveSolFreq=1;
pde.soltime = pde.saveSolFreq:pde.saveSolFreq:length(pde.dt); % steps at which solution are collected
pde.timestepOffset = 130000;
% Solver params
pde.linearsolveriter = 40;
pde.GMRESrestart = 39;
pde.linearsolvertol = 1e-5;
pde.matvectol = 1e-7;
pde.NLiter = 5;
pde.RBdim = 5;
pde.NLtol = 1e-11;

%% Mutation information
% Mutation configuration 
pde.mutationflag = 1;
pde.mutationopts = {};
pde.mutationopts{1} = "air_5";
pde.mutationopts{2} = "ChemNonEq1T";
pde.mutationopts{3} = "RRHO";
pde.mutationpath = "/Users/rloekvh/Mutationpp";

% Inflow BC: specify temperature, velocity, and mass fractions
T_inflow = 6.539171319799890e+03;
Xi_inflow = [0, 0, 0, 0.76708248854242, 0.23291751145758];
u_inflow = 5.415214474186620e+02;
% Outflow BC: specify pressure
% p_outflow = 16894.1;
p_outflow = 1.689038973365919e+04;

% Dimensional "postshock" state; used as initial conditions for inflow 
% rho_post = [0 0 0 0.00655326 0.00199076];
% rhou_post = 4.63086;
% rhoE_post = 52267.1;

% rho_post = [0.000000000000000e+00 0.000000000000000e+00 0.000000000000000e+00 6.553264087414407e-03 1.990756887050269e-03 ];
% rho_post = [0.000000000000000   0.000000000000000   0.000000000000000   0.006553622130237   0.001989946818508];
% rho_post = [6.843577361516948e-81 3.586042485049209e-42 2.417482657189227e-18 6.875039878366979e-03 2.087542348520588e-03 ];
rho_post = [6.523630530137378e-81 3.418389973844208e-42 2.304461955409462e-18 6.553622130237057e-03 1.989946818507937e-03 ];
% rhou_post = 4.630859368159855e+00;
% rhou_post = 4.626525823245757e+00;
rhou_post = 4.626525823245526e+00;
% rhoE_post = 5.226378663945410e+04;
rhoE_post = 5.226378663945407e+04;

% Scaling parameters
rho_scale = sum(rho_post);
u_scale = u_inflow;
rhoe_scale = rho_scale * u_scale * u_scale;

% Dimensional equilibrium state; used as initial conditions for outflow 
% rho_equil = [2.0979e-05 0.00259561 0.000326699 0.00941716 0.000142304];
% rhou_equil = 4.62653;
% rhoE_equil = 83151.3;
rho_equil = [2.097903177177863e-05 2.595607445498073e-03 3.266988117157821e-04 9.417156962182115e-03 1.423035901034682e-04 ];
rhou_equil = 4.626525823245757e+00;
rhoE_equil = 8.315126885726860e+04;

pde.physicsparam = [rho_post(:)', rhou_post, rhoE_post, rho_equil(:)', rhou_equil, rhoE_equil,   Xi_inflow,    u_inflow, T_inflow, p_outflow];
%                      1:N          N+1         N+2       N+3:2*N+2     2*N+3       2*N+4       2*N+5:3*N+4     3*N+5     3*N+6      3*N+7
pde.externalparam = [rho_scale, u_scale, rhoe_scale, NaN]; 
% Inf is a placeholder variable that's replaced with Mutation outputs
%% Stabilization and mesh
pde.tau = [5,5,5,5,5,5,5];
nDiv = 1024;
[mesh.p,mesh.t] = linemesh(nDiv-1);
a = 0; b = 1.5; 
% B=10;
% dx = 5e-6;
% b1 = 1e-3;
% p1 = 0:dx:(b1-dx);
% n1 = length(p1);
% mesh.p = [p1, linspace(b1, 1.5, nDiv-n1)];
mesh.p = a + (b-a)*mesh.p;
mesh.p = loginc(mesh.p, 6);
% mesh.p = [mesh.p, 1.6:0.1:B];
% m = length(mesh.p);
% mesh.t = [1:1:(m-1); 2:1:m];
disp(mesh.p(2))
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(1,:)-a)<1e-16, @(p) abs(p(1,:) - b)<1e-8}; %TODO: double check boundaries
mesh.boundarycondition = [1;2]; % Set boundary condition for each boundary

%% Call Exasim
% search compilers and set options
pde = setcompilers(pde);       

% generate input files and store them in datain folder
[pde,mesh,master,dmd] = preprocessing(pde,mesh);

% generate source codes and store them in app folder
% gencode(pde);
% cd("app");
% eval('!./fixargsMPP_Mac.sh');
% eval('!cp opuAppMPP.cpp opuApp.cpp');
% eval('!cp opuUbouSubsonicMPP.cpp opuUbou.cpp');
% % eval('!cp opuFbouSubsonicMPP.cpp opuFbou.cpp');
% eval('!cp opuFbouSubsconicVectorTau.cpp opuFbou.cpp');
% eval('!cp opuFlux_MPP.cpp opuFlux.cpp');
% eval('!cp opuSource_MPP.cpp opuSource.cpp');
% eval('!cp opuOutput_MPP.cpp opuOutput.cpp');
% cd("..");
% 
% % compile source codes to build an executable file and store it in app folder
% compilerstr = compilecode(pde);
% 
% % run executable file to compute solution and store it in dataout folder
% runstr = runcode(pde);   
