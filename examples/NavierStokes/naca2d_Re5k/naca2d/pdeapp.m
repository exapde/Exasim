% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim"); ii = ii(end);
run(cdir(1:(ii+5)) + "/install/setpath.m");

porder = 1;                     % polynomial degree
tau = 3;                      % stabilization parameter
gam = 1.4;                      % gas constant
Minf = 0.2;                     % freestream mach number
alpha = -5*pi/180;               % angle of attack
rinf = 1.0;                     % freestream density
ruinf = cos(alpha);             % freestream horizontal velocity
rvinf = sin(alpha);             % freestream vertical velocity
pinf = 1/(gam*Minf^2);          % freestream pressure
rEinf = 0.5+pinf/(gam-1);       % freestream energy
Re = 700;                       % Reynolds number 
Pr = 0.72;                      % Prandtl number 
ui = [ 1, cos(alpha), sin(alpha), 0.5+pinf/(gam-1)];

% initialize pde structure and mesh structure
[pde,~] = initializeexasim();
pde.buildpath=string(pwd()); 

pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors
pde.hybrid = 1;
pde.debugmode = 0;
pde.porder = porder;
pde.pgauss = 2*porder;

pde.physicsparam = [gam Re Pr Minf rinf ruinf rvinf rEinf];
pde.tau = tau;              % DG stabilization parameter
pde.GMRESrestart = 80;
pde.linearsolvertol = 1e-4; % GMRES tolerance
pde.linearsolveriter = 80;
pde.ppdegree = 30;
pde.RBdim = 0;
pde.neb = 512;

%% naca mesh
% mesh = mkmesh_naca0012(porder,1,10);
mesh = mkmesh_naca0012(porder,1,2);
pde.gencode=0;
pde.runmode=5;
pde.dt = 1e-4*ones(1,100);   % time step sizes
pde.soltime = [1]; % steps at which solution are collected
pde.torder = 1;          % time-stepping order of accuracy
pde.nstage = 1;          % time-stepping number of stages
pde.visdt = pde.dt(1);         % visualization timestep size
pde.saveSolFreq = 1;          % solution is saved every 10 time steps

% mesh size
[~,mesh,master,dmd] = preprocessing(pde,mesh);
mesh.dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,porder);    
mesh.porder = porder;

pde.pgauss = 2*(pde.porder);
pde.nd = 2;
pde.elemtype = 1;
master = Master(pde);
[~, ~, jac] = volgeom(master.shapent,permute(mesh.dgnodes,[1 3 2]));
hsz = reshape(sqrt(jac),[],1,size(mesh.dgnodes,3));
[~,cgelcon,rowent2elem,colent2elem,~] = mkcgent2dgent(mesh.dgnodes,1e-8);
hh = dg2cg2(max(hsz,0e-5), cgelcon, colent2elem, rowent2elem);
hh = dg2cg2(hh, cgelcon, colent2elem, rowent2elem);
mesh.vdg(:,1,:) = hh;

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh] = exasim(pde,mesh);
% [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% UDG = initu(mesh,{ui(1),ui(2),ui(3),ui(4),0,0,0,0,0,0,0,0});


%% plot solution
pde.dt = [1];
sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
figure(1); clf; scaplot(mesh,eulereval(sol(:,1:4,:),'M',gam),[]); axis off; axis equal; axis tight;
% %%
% pde.debugmode=0;
% pde.denseblock =0;
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
% % compareexasim(master, mesh1, pde);
% 
% 
% %%
% M_mat = eulereval(UDG(:,1:4,:),'M',gam);
% M_cpp = eulereval(sol(:,1:4,:),'M',gam);
% figure(1); clf; scaplot(mesh, M_mat-M_cpp);