% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
ii = ii(end);
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde structure and mesh structure
[pde,~] = initializeexasim();
pde.buildpath = string(pwd()) + "/ns";

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel_axialns";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
pde.platform = "cpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 4;              % number of MPI processors
pde.hybrid = 1;
pde.porder = 3;          % polynomial degree
pde.gencode = 1;
pde.nd = 2;

[mesh, ~] = mkmesh_orion(pde.porder);
mesh.p(2,:) = mesh.p(2,:) + 1e-3;
mesh.dgnodes(:,2,:) = mesh.dgnodes(:,2,:) + 1e-3;
L = max(mesh.p(1,:));
mesh.boundaryexpr = {@(p) abs(p(2,:)-1e-3)<1e-6, @(p) p(1,:)> L-1e-2, @(p) ((p(1,:) < -1e-3) | (p(2,:) > 2.6)), @(p) abs(p(1,:))< 20 + 1e-6};
% axis symmetric, supersonic outflow, supersonic inflow,  iso-thermal wall
mesh.boundarycondition = [4, 2, 1, 3]; % Set boundary condition for each boundary

gam = 1.4;                      % specific heat ratio
Re = 8.0e4;                       % Reynolds number
Pr = 0.71;                      % Prandtl number    
Minf = 17;                     % Mach number
Tref  = 200;
Twall = 1464;
pinf = 1/(gam*Minf^2);
Tinf = pinf/(gam-1);
alpha = 0;                % angle of attack
rinf = 1.0;                     % freestream density
ruinf = cos(alpha);             % freestream horizontal velocity
rvinf = sin(alpha);             % freestream vertical velocity
pinf = 1/(gam*Minf^2);          % freestream pressure
rEinf = 0.5+pinf/(gam-1);       % freestream energy

pde.physicsparam = [gam Re Pr Minf rinf ruinf rvinf rEinf Tinf Tref Twall];
pde.tau = 1.0;                  % DG stabilization parameter
pde.GMRESrestart = 100;         %try 50
pde.linearsolvertol = 1e-8; % GMRES tolerance
pde.linearsolveriter = 100; %try 100
pde.preconditioner = 1;
pde.GMRESortho = 1;
pde.RBdim = 0;
pde.ppdegree = 20;
pde.NLtol = 1e-6;              % Newton tolerance
pde.NLiter = 10;                % Newton iterations
pde.matvectol=1e-6;             % tolerance for matrix-vector multiplication

master = Master(pde);

% initial artificial viscosity
mesh.f = facenumbering(mesh.p,mesh.t,pde.elemtype,mesh.boundaryexpr,mesh.periodicexpr);
dist = meshdist3(mesh.f,mesh.dgnodes,master.perm,[4]); % distance to the wall
mesh.vdg = zeros(size(mesh.dgnodes,1),1,size(mesh.dgnodes,3));
mesh.vdg(:,1,:) = 0.0012*tanh(dist);
mesh.vdg(:,2,:) = 0;

mesh.porder = pde.porder;
mesh.xpe = master.xpe;
mesh.telem = master.telem;
figure(1); clf; scaplot(mesh,mesh.vdg(:,1,:),[],2); axis on; axis equal; axis tight;

load nsinit.mat 
mesh.udg = udg;
%[sol,pde,mesh,master,dmd] = exasim(pde,mesh);

[pde,mesh,master,dmd] = preprocessing(pde,mesh);
kkgencode(pde);
% compilerstr = cmakecompile(pde); % use cmake to compile C++ source codes 
% runstr = runcode(pde, 1); % run C++ code

return;



