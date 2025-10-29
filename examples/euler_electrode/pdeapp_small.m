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
pde.platform = "cpu";
pde.mpiprocs = 1;              % number of MPI processors

% Discretization parameters
pde.porder = 2;          % polynomial degree

% Solver/prec params
pde.GMRESortho = 0;   % 0 is MGS, 1 is CGS. numerical instability of CGS only kicks in around 800 and above GMRES vectors.
% pde.GMRESortho = 1;
% pde.ppdegree = 5; % only use if > 1000 gmres iterations
pde.ppdegree = 1; % only use if > 1000 gmres iterations -- polynomial preconditioner isn't that good. if you have a good preconditioner you don't really need it.
pde.linearsolvertol = 1e-4; % GMRES tolerance  -- don't need that much accuracy at the beginning of the newton iteration -- 10e-3 or 10e-4
pde.NLiter = 5;
pde.RBdim = 0;              % reduced basis dimension for preconditioner -- 5 is good
pde.linearsolveriter = 400; % tdep solver is more efficient -> fewer GMRES iterations. was 100. 200-500 for steady problems. no restarts if you can handle it
pde.GMRESrestart = 150; % larger gmres for steady problems. was 20 -- make this as high as possible, memory limited
% pde.linearsolveriter = 100; % tdep solver is more efficient -> fewer GMRES iterations. was 100
% pde.GMRESrestart = 40; % larger gmres for steady problems. was 20
pde.NLtol = 1e-7;
pde.precMatrixType = 2; % 2 is good for tdep problems, can be [0,1,2]. 0=nothing, 1=use mass matrix, 2=inv(mass matrix), good for tdep problems
pde.preconditioner = 1; % Additive Schwarz

% Timestepping
pde.torder = 2;          % time-stepping order of accuracy
pde.nstage = 2;          % time-stepping number of stages
pde.dt = 5e-3*ones(1,100);   % time step sizes
pde.soltime = 1:length(pde.dt); % steps at which solution are collected
pde.visdt = 0.05; % visualization timestep size

% Pysics param
pde.physicsparam = [1.4];     % Physics param loaded in a separate script
pde.tau = 1;
% pde.gencode = 0;

% Mesh
[p,t] = gmshcall_sam("./streamer2_39k.msh", 2, 0);

% Normalization
xmax = max(p(1,:));
p = p*10;       % Characteristic length scale is .1 mm

p = [p(2,:); p(1,:)];
mesh.p = p;
mesh.t = t;
mesh.boundaryexpr = {@(p) abs(p(2,:))<1e-3, @(p) and(p(2,:)<8, p(1,:)<0) , @(p) and(p(2,:)<8, p(1,:)>0), @(p) abs(p(2,:)-7.94)>0}; %@(p) p(2,:)>7.5};
mesh.boundarycondition = [1;2;3;4]; % Set boundary condition for each boundary
mesh.f = facenumbering(mesh.p,mesh.t,0,mesh.boundaryexpr,[]);
mesh.porder = pde.porder;
mesh.dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);

% generate input files and store them in datain folder
[pde,mesh,master,dmd] = preprocessing(pde,mesh);

% save("dmd.mat", "dmd");
% disp("Saved dmd, exiting")
% return;
mesh.plocal = mesh.xpe;


%%%%%%%%%%%% Initial condition: specify rho0 and T0
rho0 = 1.225;
T0 = 298;   % K
cv = 718;   % Perfect gas: constant specific heats
Rsp = 287;      
p0 = rho0*Rsp*T0;    % For thermodynamic consistency, we have to use the equation of state to compmute T0 from rho0 and p0; we can't specify it separately.
rE0 = rho0*cv*T0/(p0/rho0);   % Nondimensional initial rho*E is rho_0*cv(T_0)*T_0/(p_ref/rho_ref)
initu_func_set = {rho0;0;0;rE0};
UDG0 = initu_euler(mesh,initu_func_set,pde.physicsparam);

% Print out inital condition and boundary values for a verification
fig = figure(1); clf; scaplot(mesh,UDG0(:,1,:),[],1); axis on; axis equal; axis tight;
saveas(fig, "u0/U1.png")
fig = figure(1); clf; scaplot(mesh,UDG0(:,2,:),[],1); axis on; axis equal; axis tight;
saveas(fig, "u0/U2.png")
fig = figure(1); clf; scaplot(mesh,UDG0(:,3,:),[],1); axis on; axis equal; axis tight;
saveas(fig, "u0/U3.png")
fig = figure(1); clf; scaplot(mesh,UDG0(:,4,:),[],1); axis on; axis equal; axis tight;
saveas(fig, "u0/U4.png")

fig = figure(1); clf; boundaryplot(mesh,1);
saveas(fig, "u0/B1.png")
fig = figure(1); clf; boundaryplot(mesh,2);
saveas(fig, "u0/B2.png")
fig = figure(1); clf; boundaryplot(mesh,3);
saveas(fig, "u0/B3.png")
fig = figure(1); clf; boundaryplot(mesh,4);
saveas(fig, "u0/B4.png")

%%%%%%%%%%%%

mesh.udg = UDG0;            % Load initial condition

pde.saveSolFreq = 1;

% search compilers and set options
pde = setcompilers(pde);       

% return;    
[sol,pde,mesh] = exasim(pde,mesh);
