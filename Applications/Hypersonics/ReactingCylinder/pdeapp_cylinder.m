% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

pde.cpucompiler="/opt/homebrew/Cellar/llvm@12/12.0.1_1/bin/clang++";
% pde.metis = "/opt/homebrew/Cellar/metis/5.1.0/bin/mpmetis";
% pde.mpicompiler = "/usr/local/Cellar/open-mpi/4.1.2/bin/mpicxx ";
% pde.cpuflags = pde.cpuflags + " -arch x86_64";
% pde.cpulibflags = "-arch x86_64";
% pde.cpuappflags = "-arch x86_64";
% pde.cpuflags = pde.cpuflags + " -arch x86_64";

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel_2d";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
pde.mpiprocs = 1;              % number of MPI processors

%% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 2;          % polynomial degree
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
dt = 1.0000000000000e-5;
pde.dt = [dt*ones(1,5*nt)]; %, 10*dt*ones(1, nt), 100*dt*ones(1,8*nt)];
% tfinal = dt*nt + 10*dt*nt + 100*dt*8*nt
% pde.dt = [dt*ones(1,nt)];
pde.saveSolFreq = ceil(nt/100);
% pde.saveSolFreq = 100;
% pde.soltime = pde.saveSolFreq:pde.saveSolFreq:length(pde.dt); % steps at which solution are collected
pde.timestepOffset = 600;

% Solver params
pde.linearsolveriter = 80;
pde.GMRESrestart = 19;
pde.linearsolvertol = 1e-5;
pde.matvectol = 1e-7;
pde.NLiter = 4;
pde.RBdim = 40;
pde.NLtol = 1e-9;

pde.AV = 1;      
% pde.nce = 12;

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
ndim = 2;
pde.tau = [5, 5, 5, 5, 5, 5, 5, 5];

% mesh size
rbody = 0.045;
rdomain = 0.12;
[mesh.p,mesh.t,mesh.dgnodes] = mkmesh_circincirc_Ma17b(pde.porder,61,61,rbody,rdomain,1);
mesh.boundaryexpr = {@(p) sqrt(p(1,:).^2+p(2,:).^2)<rbody+1e-6, @(p) p(1,:)>-1e-7, @(p) abs(p(1,:))<20};
% adiabatic wall, supersonic outflow, supersonic inflow
mesh.boundarycondition = [3;2;1]; 

pde.pgauss = 2*(pde.porder);
pde.nd = ndim;
pde.elemtype = 1;
master = Master(pde);
[~, ~, jac] = volgeom(master.shapent,permute(mesh.dgnodes,[1 3 2]));
hsz = reshape(sqrt(jac),[],1,size(mesh.dgnodes,3));
[~,cgelcon,rowent2elem,colent2elem,~] = mkcgent2dgent(mesh.dgnodes,1e-8);
hh = dg2cg2(max(hsz,0e-5), cgelcon, colent2elem, rowent2elem);
hh = dg2cg2(hh, cgelcon, colent2elem, rowent2elem);

% distance to the wall
mesh.f = facenumbering(mesh.p,mesh.t,pde.elemtype,mesh.boundaryexpr,mesh.periodicexpr);
f = mkf(mesh.t,mesh.f,2);
dist = meshdist(f,mesh.dgnodes,master.perm,[1]); % distance to the wall

mesh.vdg = zeros(size(mesh.dgnodes,1),2,size(mesh.dgnodes,3));
mesh.vdg(:,2,:) = hh.*tanh(1000*dist);

%% Initial conditions and parameters
% TL=9000

% Dimensional inflow quantities
rho_inf = [1.0055500000000001e-09 0.00035318161606 1.5872374699999997e-05 0.00116691299088 1.103201281e-05];
rhou_inf = 9.213932;
rhov_inf = 0.0;
rhoE_inf = 33563.20282790763;

% Reference quantities for nondimensionalizing
rho_ref = 0.001547;
% u_ref = 665.5369622053455;
u_ref = 5956;
% rhoE_ref = 685.2273261511706;
rhoE_ref = rho_ref * u_ref^2;
T_ref = 374.6764014579334;
mu_ref = 4.415107627874723e-05;
kappa_ref = 0.07663957108841235;

% Nondimensional quantities 
Re = 23319.605484391683;
Pr = 0.681045690355942;

% Av params
kb = 1.5;
sb0 = 0.001;
sbmax = 2.5;

% Load into Exasim data structures
pde.physicsparam(1:nspecies) = rho_inf;
pde.physicsparam(nspecies+1) = rhou_inf;
pde.physicsparam(nspecies+2) = rhov_inf;
pde.physicsparam(nspecies+3) = rhoE_inf;
pde.physicsparam(15) = pde.porder;
pde.physicsparam(16) = kb;
pde.physicsparam(17) = sb0;
pde.physicsparam(18) = sbmax;
pde.physicsparam(19) = Re;
pde.physicsparam(20) = Pr;

pde.externalparam = zeros(3*nspecies + 3 + ndim);
pde.externalparam(1) = rho_ref; % rho_inf
pde.externalparam(2) = u_ref;    % u_inf
pde.externalparam(3) = rhoE_ref;  % rhoE_inf
pde.externalparam(4) = T_ref;
pde.externalparam(5) = mu_ref;
pde.externalparam(6) = kappa_ref;

% intial solution
ui = [rho_inf(:)'/rho_ref rhou_inf/(rho_ref*u_ref) rhov_inf/(rho_ref*u_ref) rhoE_inf/rhoE_ref];
UDG = initu(mesh,{ui(1),ui(2),ui(3),ui(4),ui(5),ui(6),ui(7),ui(8),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}); % freestream 
UDG(:,nspecies+1,:) = UDG(:,nspecies+1,:).*tanh(100*dist);
UDG(:,nspecies+2,:) = UDG(:,nspecies+2,:).*tanh(100*dist);
mesh.udg = UDG;

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
eval('!cp opuFlux_2d_MPP.cpp opuFlux.cpp'); %TODO
eval('!cp opuFbou_2d_MPP.cpp opuFbou.cpp'); %TODO
eval('!cp opuSource_2d_MPP.cpp opuSource.cpp'); %X
eval('!cp opuOutput_2d_MPP.cpp opuOutput.cpp'); %TODO
eval('!cp opuAvfield_2d_MPP.cpp opuAvfield.cpp') %X
% eval('!cp opuUbou_MPP.cpp opuUbou.cpp')
cd("..");

% compile source codes to build an executable file and store it in app folder
compilerstr = compilecode(pde);

% run executable file to compute solution and store it in dataout folder
runstr = runcode(pde);   
