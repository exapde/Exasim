% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
ii = ii(end);
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde structure and mesh structure
tmp = load('meshIDG.mat');
meshIDG = tmp.mesh;
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel_transonic_2";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
%pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.cpucompiler="/opt/homebrew/Cellar/llvm@12/12.0.1_1/bin/clang++";
% pde.enzyme = "ClangEnzyme-12.dylib";
% pde.metis = "/opt/homebrew/Cellar/metis/5.1.0/bin/mpmetis";
% pde.mpicompiler = "/usr/local/Cellar/open-mpi/4.1.2/bin/mpicxx";

% pde.cpuflags = pde.cpuflags + " -arch x86_64";
% pde.cpulibflags = "-arch x86_64";
% pde.cpuappflags = "-arch x86_64";

pde.mpiprocs = 1;              % number of MPI processors

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 2;          % polynomial degree
pde.torder = 1;          % time-stepping order of accuracy
pde.nstage = 1;          % time-stepping number of stages
pde.dt = 4*0.25*[2*0.0003*ones(1,100), 4*0.0003*ones(1,100), 4*0.0003*ones(1,200), 4*0.0003*ones(1,10*700)];
% pde.timestepOffset = 7400;
% pde.timestepOffset = 3200;
% pde.dt = [0.0003*ones(1,100),0.03*ones(1,100),0.3*ones(1,100)]
pde.saveSolFreq = 100;          % solution is saved every 10 time steps
pde.soltime = pde.saveSolFreq:pde.saveSolFreq:length(pde.dt);% steps at which solution are collected
pde.visdt = pde.dt(1);           % visualization timestep size

gam = 1.4;              % specific heat ratio
Minf = 0.8;             % Mach number
rinf = 1.0;             % freestream density
aoa = -4;
alpha = aoa*pi/180;
uinf = cos(alpha);             % freestream horizontal velocity
vinf = sin(alpha);             % freestream vertical velocity
pinf = 1/(gam*Minf^2);  % freestream pressure
rEinf = 0.5+pinf/(gam-1); % freestream energy

avb = 1.5;                     % bulk viscosity parameter
sb0 = 0.02;                     % cutoff  dilatation
sb1 = 2.5;    

pde.physicsparam = [gam Minf rinf uinf vinf rEinf, avb, pde.porder, sb0, sb1];

pde.tau = 2.0;          % DG stabilization parameter

pde.GMRESrestart=30;  % number of GMRES restarts
pde.linearsolvertol=0.01; % GMRES tolerance
pde.linearsolveriter=31;  % number of GMRES iterations
pde.precMatrixType=2; % preconditioning type
pde.NLtol = 1e-7;  % Newton tolerance
% pde.timestepOffset = 50000;
pde.NLiter = 2;   % number of Newton iterations
% pde.ppdegree = 15;
pde.AV = 1;

% read a grid from a file
% [mesh.p,mesh.t] = readmesh('grid_ns.bin',0); pde.elemtype = 1;

mesh.p = meshIDG.p';
mesh.t = meshIDG.t'; pde.elemtype = 0;
% [mesh.p,mesh.t] = uniref(mesh.p, mesh.t, 1);
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) sqrt((p(1,:)-.5).^2+p(2,:).^2)<3, @(p) abs(p(1,:))<20};
% mesh.boundarexpr = {@(p) sqrt(p(1,:).^2 + p(2,:).^2)<2,@(p) sqrt(p(1,:).^2 + p(2,:).^2)<20};  
mesh.boundarycondition = [2;1];
% expressions for curved boundaries
mesh.curvedboundary = [1 0];
mesh.curvedboundaryexpr = {@(p) p(2,:).^2-(5*0.01*12*(0.29690*sqrt(abs(p(1,:)))-0.12600*p(1,:)-0.35160*p(1,:).^2+0.28430*p(1,:).^3-0.10150*p(1,:).^4)).^2, @(p) 0};

% mesh size
[pde,mesh,master,dmd] = preprocessing(pde,mesh);
dgnodes = mesh.dgnodes;
pde.pgauss = 2*(pde.porder);
pde.nd = 2;
master = Master(pde);
[~, ~, jac] = volgeom(master.shapent,permute(mesh.dgnodes,[1 3 2]));
hsz = reshape(sqrt(jac),[],1,size(mesh.dgnodes,3));
[~,cgelcon,rowent2elem,colent2elem,~] = mkcgent2dgent(mesh.dgnodes,1e-8);
hh = dg2cg2(max(hsz,0e-5), cgelcon, colent2elem, rowent2elem);
hh = dg2cg2(hh, cgelcon, colent2elem, rowent2elem);
% hh = 0*hh + 0.01;

% distance to the wall
mesh.f = facenumbering(mesh.p,mesh.t,pde.elemtype,mesh.boundaryexpr,mesh.periodicexpr);
f = mkf(mesh.t,mesh.f,2);
dist = meshdist(f,mesh.dgnodes,master.perm,[1]); % distance to the wall

mesh.vdg = zeros(size(mesh.dgnodes,1),2,size(mesh.dgnodes,3));
mesh.vdg(:,2,:) = hh.*tanh(1e6*dist);

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh,master,dmd,compilerstr,runstr,res] = exasim(pde,mesh);
% pde = setcompilers(pde);       
% 
% % generate input files and store them in datain folder
% [pde,mesh,master,dmd] = preprocessing(pde,mesh);

% generate source codes and store them in app folder
% gencode(pde);
%%
% visualize the numerical solution of the PDE model using Paraview
pde.visscalars = {"density", 1, "energy", 4};  % list of scalar fields for visualization
pde.visvectors = {"momentum", [2, 3]}; % list of vector fields for visualization
xdg = vis(sol,pde,mesh); % visualize the numerical solution
disp("Done!");
