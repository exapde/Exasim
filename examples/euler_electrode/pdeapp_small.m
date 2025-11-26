% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";       % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel_axialns";    % name of a file defining the PDE model
% pde.buildpath = string(pwd()) + "/ns";

pde.hybrid = 1;
pde.debugmode = 0;

% Choose computing platform and set number of processors
pde.platform = "gpu";
pde.mpiprocs = 1;              % number of MPI processors

% Discretization parameters
pde.porder = 2;          % polynomial degree

% Solver/prec params
pde.GMRESortho = 0;   % 0 is MGS, 1 is CGS. numerical instability of CGS only kicks in around 800 and above GMRES vectors.
% pde.GMRESortho = 1;
% pde.ppdegree = 5; % only use if > 1000 gmres iterations
pde.ppdegree = 0; % only use if > 1000 gmres iterations -- polynomial preconditioner isn't that good. if you have a good preconditioner you don't really need it.
pde.linearsolvertol = 1e-7; % GMRES tolerance  -- don't need that much accuracy at the beginning of the newton iteration -- 10e-3 or 10e-4
pde.NLiter = 5;
pde.RBdim = 5;              % reduced basis dimension for preconditioner -- 5 is good
pde.linearsolveriter = 400; % tdep solver is more efficient -> fewer GMRES iterations. was 100. 200-500 for steady problems. no restarts if you can handle it
pde.GMRESrestart = 200; % larger gmres for steady problems. was 20 -- make this as high as possible, memory limited
% pde.linearsolveriter = 100; % tdep solver is more efficient -> fewer GMRES iterations. was 100
% pde.GMRESrestart = 40; % larger gmres for steady problems. was 20
pde.NLtol = 1e-8;
pde.precMatrixType = 2; % 2 is good for tdep problems, can be [0,1,2]. 0=nothing, 1=use mass matrix, 2=inv(mass matrix), good for tdep problems
pde.preconditioner = 1; % Additive Schwarz

% pde.gencode=0;
% Timestepping
pde.torder = 2;          % time-stepping order of accuracy
pde.nstage = 2;          % time-stepping number of stages
pde.dt = 1e-3*ones(5000,1);   % time step sizes
pde.soltime = 1:length(pde.dt); % steps at which solution are collected
pde.visdt = 0.05; % visualization timestep size

% Pysics param
pde.physicsparam = [1.4, 0e-7];     % Physics param loaded in a separate script
pde.tau = 0.01;
pde.tau = 1e-3*ones(4,1);
% pde.tau(1) = 1e-4;
% pde.tau(2) = 1e-4;
% pde.tau(3) = 1e-4;
% pde.tau(4) = 1e-4;
% pde.tau(4) = 1e-4;
% pde.tau(4) = 1e-4;
pde.gencode = 1;

% Mesh
[p,t] = gmshcall_sam("./euler_mesh_390.msh", 2, 0);
% [p,t] = gmshcall_sam("./streamer2_39k.msh", 2, 0);

% Normalization
xmax = max(p(1,:));
p = p*10;       % Characteristic length scale is .1 mm

% p = [p(2,:); p(1,:)];     % Do not transpose the 390K element mesh
mesh.p = p;
mesh.t = t;
mesh.boundaryexpr = {@(p) abs(p(2,:))<1e-3, @(p) and(p(2,:)<8, p(1,:)<0) , @(p) and(p(2,:)<8, p(1,:)>0), @(p) abs(p(2,:)-7.94)>0}; %@(p) p(2,:)>7.5};
mesh.boundarycondition = [2;2;3;4]; % Set boundary condition for each boundary
mesh.f = facenumbering(mesh.p,mesh.t,0,mesh.boundaryexpr,[]);
mesh.porder = pde.porder;
mesh.dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);

% generate input files and store them in datain folder
% [pde,mesh,master,dmd] = preprocessing(pde,mesh);

% save("dmd.mat", "dmd");
% disp("Saved dmd, exiting")
% return;
%mesh.plocal = mesh.xpe;
[mesh.xpe,mesh.telem,xpf,tface,master.perm] = masternodes(pde.porder,2,0);
mesh.plocal = mesh.xpe;

%%%%%%%%%%%% Initial (reference) condition: specify two of rho, T, p, or E. In this simulation I have chosen to specify rho and T.
rho0 = 1.225;
T0 = 298;   % K
cv = 718;   % Perfect gas: constant specific heats
Rsp = 287;
p0 = rho0*Rsp*T0;
rE0 = rho0*cv*T0/p0;   % Nondimensional initial rho*E = 2.5
% initu_func_set = {rho0;0;0;rE0;};
initu_func_set = {rho0;0;0;@initu_func_rE;};
UDG0 = initu_euler(mesh,initu_func_set,pde.physicsparam);

% % Print out inital condition and boundary values for a verification
% fig = figure(1); clf; scaplot(mesh,UDG0(:,1,:),[],1); axis on; axis equal; axis tight;
% saveas(fig, "u0/U1.png")
% fig = figure(2); clf; scaplot(mesh,UDG0(:,2,:),[],1); axis on; axis equal; axis tight;
% saveas(fig, "u0/U2.png")
% fig = figure(3); clf; scaplot(mesh,UDG0(:,3,:),[],1); axis on; axis equal; axis tight;
% saveas(fig, "u0/U3.png")
% fig = figure(4); clf; scaplot(mesh,UDG0(:,4,:),[],1); axis on; axis equal; axis tight;
% saveas(fig, "u0/U4.png")

% fig = figure(5); clf; boundaryplot(mesh,1);
% saveas(fig, "u0/B1.png")
% fig = figure(7); clf; boundaryplot(mesh,2);
% saveas(fig, "u0/B2.png")
% fig = figure(7); clf; boundaryplot(mesh,3);
% saveas(fig, "u0/B3.png")
% fig = figure(8); clf; boundaryplot(mesh,4);
% saveas(fig, "u0/B4.png")

%%%%%%%%%%%%

%UDG0 = sol0(:,:,:,50);
mesh.udg = UDG0;            % Load initial condition

pde.saveSolFreq = 10;

% search compilers and set options
pde = setcompilers(pde);

mesh.f = facenumbering(mesh.p,mesh.t,0,mesh.boundaryexpr,mesh.periodicexpr);
dist = meshdist3(mesh.f,mesh.dgnodes,master.perm,[1 2 3]); % distance to the wall
mesh.vdg = zeros(size(mesh.dgnodes,1),1,size(mesh.dgnodes,3));
% mesh.vdg(:,1,:) = 0; %0e-6*tanh(dist);
AV = 1e-4;
mesh.vdg(:,1,:) = AV;
disp("Assigning viscosity")
disp(AV)

% return;
[sol,pde,mesh] = exasim(pde,mesh);

%[k,sol] = readsoldmd("/Users/cuongnguyen/Documents/GitHub/Exasim/build/dataout/outudg", dmd, 1, 45);
%figure(4); clf; scaplot(mesh,sol(:,1,:),[],1); axis on; axis equal; axis tight;

% figure(2); clf; scaplot(mesh,sol(:,2,:,end)./sol(:,1,:,end),[],1);
% axis on; axis equal; axis tight;
% figure(3); clf; scaplot(mesh,sol(:,3,:,end)./sol(:,1,:,end),[],1);
% axis on; axis equal; axis tight;
% figure(4); clf; scaplot(mesh,sol(:,1,:,end),[],1);
% axis on; axis equal; axis tight;
