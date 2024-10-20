% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

porder = 4;                     % polynomial degree
tau = 3.0;                      % stabilization parameter
gam = 1.4;                      % gas constant
Minf = 0.2;                     % freestream mach number
alpha = 0*pi/180;               % angle of attack
rinf = 1.0;                     % freestream density
ruinf = cos(alpha);             % freestream horizontal velocity
rvinf = sin(alpha);             % freestream vertical velocity
pinf = 1/(gam*Minf^2);          % freestream pressure
rEinf = 0.5+pinf/(gam-1);       % freestream energy
Re = 100;                       % Reynolds number 
Pr = 0.72;                      % Prandtl number 
ui = [ 1, cos(alpha), sin(alpha), 0.5+pinf/(gam-1)];

% initialize pde structure and mesh structure
[pde,~] = initializeexasim();

pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
%pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors
pde.hybrid = 1;
pde.debugmode = 0;
pde.porder = porder;
pde.pgauss = 2*porder;

pde.physicsparam = [gam Re Pr Minf rinf ruinf rvinf rEinf];
pde.tau = tau;              % DG stabilization parameter
pde.GMRESrestart = 400;
pde.linearsolvertol = 1e-7; % GMRES tolerance
pde.linearsolveriter = 400;
pde.ppdegree = 0;
pde.RBdim = 0;
pde.neb = 512;

% naca mesh
mesh = mkmesh_naca0012(porder,1,1);

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh] = exasim(pde,mesh);

% plot solution
mesh.porder = porder;
mesh.dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,porder);    
figure(1); clf; scaplot(mesh,eulereval(sol(:,1:4,:),'M',gam),[],2); axis off; axis equal; axis tight;

% pde.denseblock = 1;
% pde.source = 'source';
% pde.flux = 'flux';
% pde.fbou = 'fbou';
% pde.fhat = 'fhat';
% pde.arg = {gam,0.0,Re,Pr,Minf,tau};
% pde.bcm  = [2,1];  % 2: Wall, 1: Far-field
% pde.bcs  = [ui;ui];
% pde.fc_q = 1;
% pde.fc_u = 1;
% pde.time = 0;
% pde.tdep = 0;
% pde.nd = 2;
% mesh1 = hdgmesh(mesh, porder);
% master = Master(pde);
% UDG0 = initu(mesh1,{ui(1),ui(2),ui(3),ui(4); 0,0,0,0; 0,0,0,0});
% UH0 = getuhat(UDG0, mesh1.f2t, master.perm, 4);
% [UDG,UH] = hdgsolve(master,mesh1,pde,UDG0,UH0,[]);
% 
% compareexasim(master, mesh1, pde);


