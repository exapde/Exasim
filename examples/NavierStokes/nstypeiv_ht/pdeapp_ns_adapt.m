cdir = pwd(); ii = strfind(cdir, "Exasim");
ii = ii(end);
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde_ns structure and mesh structure
[pde_ns,~] = initializeexasim();
pde_ns.buildpath = string(pwd()) + "/ns";

% Define a pde_ns model: governing equations, initial solutions, and boundary conditions
pde_ns.model = "ModelD";          % ModelC, ModelD, ModelW
pde_ns.modelfile = "pdemodel_ns";    % name of a file defining the pde_ns model

% Choose computing platform and set number of processors
pde_ns.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde_ns.mpiprocs = 1;              % number of MPI processors
pde_ns.hybrid = 1;
pde_ns.porder = 4;          % polynomial degree
pde_ns.gencode = 1;

porder = 4;
hybrid = 'hdg'; 

gam = 1.4;
Minf = 8.03;                  % Infinity conditions
Re = 1.94e5;
Pr = 0.71;
tau = 4;
kk = 30;
Tref  = 122.11;
Twall = 294.44;
nd    = 2;

gam1 = gam-1.0;
Tinf = 1/(gam*gam1*Minf^2);
TisoW = (Twall/Tref) * Tinf;

epslm = 0.0;
tau = 4;

pinf = 1/(gam*Minf^2);
ui = [ 1.0, 1.0, 0.0, 0.5+pinf/(gam-1)];

pde_ns.hybrid = hybrid;
pde_ns.localsolve = 1;
pde_ns.arg = {gam,Minf,epslm,Re,Pr,Tref,Twall,kk,tau};
pde_ns.physicsparam = [gam Re Pr Minf 0 0 0 0 Tinf Tref Twall 0 0 0 pde_ns.porder 0 0];

pde_ns.bcm  = [8, 5, 9];  
pde_ns.bcs  = [ui; ui; ui; ui]; 

pde_ns.tdep = false;
pde_ns.wave = false;
pde_ns.alag = false;
pde_ns.flg_q = 1;
pde_ns.flg_p = 0;
pde_ns.flg_g = 0;

pde_ns.ndim = 2;
pde_ns.nch  = 2+pde_ns.ndim;                % Number of componets of UH
pde_ns.nc   = pde_ns.nch*3;                 % Number of componeents of UDG
pde_ns.ncu  = pde_ns.nch;                   % Number of components of U

pde_ns.time = [];

pde_ns.fc_q = 1;
pde_ns.fc_u = 0;
pde_ns.bcs  = [ui; ui; ui; ui];
pde_ns.bcd  = [1,1,1,1];  
pde_ns.bcv  = [0; 0; 0; 0];

%TODO: need to get the mesh gen in there. 
mesh_ns = mkmesh_square2(61,61,porder,1,1,1,1,1);
mesh_ns = mkmesh_halfcircle(mesh_ns, 1, 3, 4.7, pi/2, 3*pi/2);

mesh_ns.porder = porder;
mesh_ns.boundaryexpr = {@(p) sqrt(p(1,:).^2+p(2,:).^2)<1+1e-6, @(p) p(1,:)>-1e-7, @(p) abs(p(1,:))<20};
mesh_ns.bndexpr = {'all(sqrt(p(:,1).^2+p(:,2).^2)<1+1e-6)'    'all(p(:,1)>-1e-7)'    'all(abs(p(:,1))<20)'};
mesh_ns.periodicexpr = {};
% clf; meshplot(mesh_ns);
mesh_ns.boundarycondition = [8;5;9]; 


%  TODO: need dist
% mesh_ns.dist = tanh(meshdist(mesh_ns,2)*10);
dist = load('dist.mat','dist');
mesh_ns.dist = dist;

[pde_ns,mesh,master,dmd] = preprocessing(pde_ns,mesh_ns);

% kappa = kappalist(jiter);
% tau = 4;

mat_res = load('result.mat');


pde_ns.source = 'source_ns';
pde_ns.flux = 'flux3';
pde_ns.fbou = 'fbou_ns';
pde_ns.fhat = 'fhat3';

% pde_ns.arg = {gam,0.0,Re,Pr,Minf,tau};
% pde_ns.bcm  = [2,1];  % 2: Wall, 1: Far-field
% pde_ns.bcs  = [ui;ui];
pde_ns.fc_q = 1;
pde_ns.fc_u = 1;
pde_ns.time = 0;
pde_ns.tdep = 0;

pde_ns.nd = 2;
mesh1 = hdgmesh(mesh_ns, porder);
master = Master(pde_ns);


%%% Get working on adapted mesh
p2 = dgnodes_to_p(mat_res.mesh2.dgnodes(:,1:2,:), mesh1, master);
t2 = mat_res.mesh2.t;

mesh2 = mkmesh(p2,t2,porder,mesh1.bndexpr,1,1);
mesh2.p = mesh2.p';
mesh2.t = mesh2.t';
mesh2.dgnodes = mat_res.mesh2.dgnodes(:,1:2,:);

mesh2.porder = porder;
mesh2.boundaryexpr = {@(p) sqrt(p(1,:).^2+p(2,:).^2)<1+1e-4, @(p) p(1,:)>-1e-3, @(p) abs(p(1,:))<20};
mesh2.bndexpr = {'all(sqrt(p(:,1).^2+p(:,2).^2)<1+1e-4)'    'all(p(:,1)>-1e-3)'    'all(abs(p(:,1))<20)'};
mesh2.periodicexpr = {};

mesh2.boundarycondition = [8;5;9]; 

[pde_ns,mesh2,master,dmd] = preprocessing(pde_ns,mesh2);

mesh_adapt = hdgmesh(mesh2, porder);
UDG0 = mat_res.UDG2{5};
UH0  = getuhat(UDG0, mesh_adapt.f2t, master.perm, 4);
a          = mat_res.ACG2{5};
mesh_adapt.dgnodes(:,3,:) = a;

[UDG,UH] = hdgsolve(master,mesh_adapt,pde_ns,UDG0,UH0,[]);
