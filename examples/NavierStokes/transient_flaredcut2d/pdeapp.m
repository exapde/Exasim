% This example requires (mesh, master, sol) from /examples/NavierStokes/flaredplate2d
% First, you need to run /examples/NavierStokes/flaredplate2d/pdeapp.m
% Next, uncomment and run the below lines to obtain sol1 and av1

porder = 3;
mesh1 = mkmesh_flatcase2d(porder);
sol1 = fieldatdgnodes(mesh, master, sol, mesh1.dgnodes);
av1 = fieldatdgnodes(mesh, master, mesh.vdg, mesh1.dgnodes);

% Finally, you run this script

% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde structure and mesh structure
pde = initializeexasim();
pde.model = "ModelD";  
pde.modelfile = "pdemodel";

% Choose computing platform and set number of processors
pde.platform = "cpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 2;             % number of MPI processors
pde.porder = 3;          % polynomial degree
pde.pgauss = 2*pde.porder;
pde.hybrid = 1;               % 0 -> LDG, 1 -> HDG
pde.debugmode = 0;
pde.nd = 2;

gam = 1.4;                      % specific heat ratio
Re = 9.2e6;                       % Reynolds number
Pr = 0.71;                      % Prandtl number    
Minf = 7;                     % Mach number
Tref  = 51.859;
Twall = 296;
pinf = 1/(gam*Minf^2);
Tinf = pinf/(gam-1);
alpha = 0;                % angle of attack
rinf = 1.0;                     % freestream density
ruinf = cos(alpha);             % freestream horizontal velocity
rvinf = sin(alpha);             % freestream vertical velocity
pinf = 1/(gam*Minf^2);          % freestream pressure
rEinf = 0.5+pinf/(gam-1);       % freestream energy

pde.physicsparam = [gam Re Pr Minf rinf ruinf rvinf rEinf Tinf Tref Twall];
pde.tau = 2.0;                  % DG stabilization parameter
pde.GMRESrestart = 250;         %try 50
pde.GMRESortho = 1;
pde.linearsolvertol = 1e-6; % GMRES tolerance
pde.linearsolveriter = 500; %try 100
pde.preconditioner = 1;
pde.RBdim = 0;
pde.ppdegree = 0;
pde.NLtol = 1e-7;              % Newton tolerance
pde.NLiter = 10;                % Newton iterations
pde.matvectol=1e-6;             % tolerance for matrix-vector multiplication


%temporal discretization 
pde.dt = 1e-4*ones(1,20000);   % time step sizes
pde.soltime = 1:length(pde.dt); % steps at which solution are collected
%pde.visdt = 0.1; % visualization timestep size
pde.nstage = 2;
pde.torder = 2;
pde.saveSolFreq = 10;


mesh = mkmesh_flatcase2d(pde.porder);
%mesh = mkmesh_flatcase2d(pde.porder, 131, 18, 71, 0.0008);
% Inflow (left boundary), Wall (bottom boundary), Outflow (right boundary), Farfield (top boundary)
%mesh.boundarycondition = [1 4 3 2]; 
mesh.boundarycondition = [1000 4 3 2]; 

stgNmode = 200;
gridLength = 1e-03;
turbLengthFactor = 10;
visc = 3.3551679839273703e-05/2;
turbIntensity = 1/100;
Ustg = 1010.38;
pde.stgdata = stghomogeneousturbulence(gridLength, turbLengthFactor, visc, turbIntensity, Ustg, stgNmode+1);
pde.stgparam = [1 0 0];
pde.stgNmode = stgNmode;

master = Master(pde);
mesh.porder = pde.porder;
mesh.xpe = master.xpe;
mesh.telem = master.telem;

mesh.udg = sol1; % -> (u, q), which is the initial guess of the solution
mesh.vdg = sol1(:,1:4,:); % -> v(1:4) which is used to specify the turbulence inflow condition 
mesh.vdg(:,5,:) = av1;  % -> v(1) which is the AV field


[pde,mesh,master,dmd] = preprocessing(pde,mesh);%for solution extraction
%[sol,pde,mesh,master,dmd] = exasim(pde,mesh);
figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],2);
figure(2); clf; scaplot(mesh, mesh.vdg(:,5,:),[],2);

x = mesh.dgnodes(:,1,:);
y = mesh.dgnodes(:,2,:);
r = sol(:,1,:);
u = sol(:,2,:)./r;
v = sol(:,3,:)./r;

xmin = min(x(:));
ind = abs(x(:)-xmin)<1e-6;
y0 = y(ind);
u0 = u(ind);
v0 = v(ind);
r0 = r(ind);

figure(3); plot(y0, u0);
figure(4); plot(y0, v0);
figure(5); plot(y0, r0);




