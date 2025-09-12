% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim"); ii = ii(end);
run(cdir(1:(ii+5)) + "/install/setpath.m");

porder = 1;                     % polynomial degree
tau = 2;                      % stabilization parameter
gam = 1.4;                      % gas constant
Minf = 0.2;                     % freestream mach number
alpha = 0*pi/180;               % angle of attack
rinf = 1.0;                     % freestream density
ruinf = cos(alpha);             % freestream horizontal velocity
rvinf = sin(alpha);             % freestream vertical velocity
pinf = 1/(gam*Minf^2);          % freestream pressure
rEinf = 0.5+pinf/(gam-1);       % freestream energy
Re = 5000;                       % Reynolds number 
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
pde.GMRESrestart = 400;
pde.linearsolvertol = 1e-3; % GMRES tolerance
pde.NLtol = 1e-8;
pde.linearsolveriter = 400;
pde.ppdegree = 30;
pde.RBdim = 0;
pde.neb = 512;

% naca mesh
% mesh = mkmesh_naca0012(porder,1,10);
mesh = mkmesh_naca0012(porder,1,10);
pde.runmode=10; % runmode 10: timesteps until steady residual is below NLtol
% pde.runmode=11; % runmode 11: timesteps until small change in monitor function, then steady solve

pde.dt = 1e-2*ones(1,500);   % time step sizes
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

% % call exasim to generate and run C++ code to solve the PDE model
% [sol,pde,mesh] = exasim(pde,mesh);
% % [pde,mesh,master,dmd] = preprocessing(pde,mesh);
UDG = initu(mesh,{ui(1),ui(2),ui(3),ui(4),0,0,0,0,0,0,0,0});
pde.gencode=1;

disp("~~~~~~~~~~~~~~~~~~~~")
disp("~~~~~~~~~~P1~~~~~~~~~")
disp("~~~~~~~~~~~~~~~~~~~~")
UDG1 = pdeapp_porder_func(pde, mesh, porder, UDG);
%%
% pde.dt = [0];
% pde.runmode = 0;
pde.gencode=0;
pde.dt = 100*pde.dt;
for porder = 2:5
    disp("~~~~~~~~~~~~~~~~~~~~")
    disp("~~~~~~~~~~P"+string(porder)+"~~~~~~~~~")
    disp("~~~~~~~~~~~~~~~~~~~~")
    [UDG1, mesh] = pdeapp_porder_func(pde, mesh, porder, UDG1, 1);
end