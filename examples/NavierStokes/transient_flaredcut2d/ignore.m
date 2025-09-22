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
Re = 9.0996e6;                       % Reynolds number
Pr = 0.71;                      % Prandtl number    
Minf = 7;                     % Mach number
Tref  = 51.85;
Twall = 296;
pinf = 1/(gam*Minf^2);
Tinf = pinf/(gam-1);
alpha = 0;                % angle of attack
rinf = 1.0;                     % freestream density
ruinf = cos(alpha);             % freestream horizontal velocity
rvinf = sin(alpha);             % freestream vertical velocity
pinf = 1/(gam*Minf^2);          % freestream pressure
rEinf = 0.5+pinf/(gam-1);       % freestream energy

pde.Cxxpreprocessing = 0; %gpu fix
pde.physicsparam = [gam Re Pr Minf rinf ruinf rvinf rEinf Tinf Tref Twall];
pde.tau = 2.0;                  % DG stabilization parameter
pde.GMRESrestart = 250;         %orig 250
pde.GMRESortho = 1;
pde.linearsolvertol = 1e-6; % GMRES tolerance
pde.linearsolveriter = 500; %try 100
pde.preconditioner = 1;
pde.RBdim = 0;
pde.ppdegree = 0;
pde.NLtol = 1e-7;              % Newton tolerance
pde.NLiter = 10;                % Newton iterations
pde.matvectol=1e-6;             % tolerance for matrix-vector multiplication

%transient parameters
pde.dt = 5e-5*ones(1,7500);   % time step sizes
pde.soltime = 1:length(pde.dt); % steps at which solution are collected
pde.visdt = 0.1; % visualization timestep size

%turbulence parameters 
%stgdata = freestreamturb;
mesh = mkmesh_flatcase2d(pde.porder);
%mesh = mkmesh_flatcase2d(pde.porder, 131, 18, 71, 0.0008);


% Synthetic Turbulent Inflow (left boundary), Wall (bottom boundary), Outflow (right boundary), Farfield (top boundary)
mesh.boundarycondition = [1000 4 3 2]; 
pde.stgdata = stgdata;
pde.stgparam = [1 0 0];
N=200;
pde.stgNmode = N;


master = Master(pde);

mesh.porder = pde.porder;
mesh.xpe = master.xpe;
mesh.telem = master.telem;

mesh.udg = sol1; % -> (u, q), which is the initial guess of the solution
mesh.vdg = av1;  % -> v(1) which is the AV field % unchanged
mesh.vdg(:,2:5,:) = sol1(:,1:4,:); % -> v(2:5) which is used to define the inflow condition f_inflow = v(2:5) - uhat; 

%[pde,mesh,master,dmd] = preprocessing(pde,mesh);%for solution extraction

[sol,pde,mesh,master,dmd] = exasim(pde,mesh);
figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],2);
figure(2); clf; scaplot(mesh, mesh.vdg(:,1,:),[],2);
