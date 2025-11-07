% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

porder = 3;                     % polynomial degree
gam = 1.4;                      % gas constant
Minf = 0.025;                   % freestream mach number
tau = 0.6/Minf;                 % stabilization parameter
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
pde.platform = "cpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors
pde.hybrid = 1;
pde.debugmode = 0;
pde.porder = porder;
pde.pgauss = 2*porder;

pde.physicsparam = [gam Re Pr Minf rinf ruinf rvinf rEinf];
pde.tau = tau;              % DG stabilization parameter
pde.GMRESrestart = 100;
pde.GMRESortho = 1;
pde.linearsolvertol = 1e-6; % GMRES tolerance
pde.linearsolveriter = 1000;
pde.preconditioner = 1;
pde.NLtol = 1e-8;
pde.ppdegree = 10;
pde.RBdim = 0;
pde.gencode = 1;

% naca mesh
mesh = mkmesh_naca0012(porder,1,2);

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh] = exasim(pde,mesh);

% plot solution
mesh.porder = porder;
mesh.dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,porder);    

figure(1); clf; scaplot(mesh,eulereval(sol(:,1:4,:),'M',gam),[],2); 
axis on; axis equal; axis tight;
axis([-0.5 2 -0.62 0.62])
set(gca,'fontsize', 16);
%exportgraphics(gca,"nacap3M0025.png",'Resolution',200);








% pde.codegenerator = "text2code";
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


