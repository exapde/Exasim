% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";       % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model
pde.hybrid = 1;
pde.debugmode = 0;

% Choose computing platform and set number of processors
pde.platform = "gpu";
pde.mpiprocs = 1;              % number of MPI processors

% Discretization parameters
pde.porder = 2;          % polynomial degree

% Solver/prec params
pde.GMRESortho = 0;
pde.ppdegree = 0;
pde.linearsolvertol = 1e-6; % GMRES tolerance
pde.NLiter = 5;
pde.RBdim = 5;              % reduced basis dimension for preconditioner
pde.linearsolveriter = 100; % tdep solver is more efficient -> fewer GMRES iterations
pde.GMRESrestart = 20; % larger gmres for steady problems
pde.NLtol = 1e-7;
pde.precMatrixType = 2; % 2 is good for tdep problems, can be [0,1,2]. 0=nothing, 1=use mass matrix, 2=inv(mass matrix), good for tdep problems

% Timestepping
pde.torder = 2;          % time-stepping order of accuracy
pde.nstage = 2;          % time-stepping number of stages
pde.dt = 5e-3*ones(1,200);   % time step sizes
pde.soltime = 1:length(pde.dt); % steps at which solution are collected
pde.visdt = 0.05; % visualization timestep size

% Pysics param
ne_star = 10^7; % this constant is used to normalize ne and ni
pde.physicsparam = phys_param();     % Physics param loaded in a separate script
pde.physicsparam(end+1) = ne_star;
pde.tau = 1;            % DG stabilization parameter

% Mesh
[p,t] = gmshcall("streamer_98k.msh", 2, 0);
% [p,t] = gmshcall("streamer_350k.msh", 2, 0);
xmax = max(p(1,:));
p = p/xmax * 125;   % Normalization
mesh.p = p;
mesh.t = t;

% Before running this line, call "mesh1 =n mkmesh_streamer_gmsh(2,"streamer_16k_fixed.msh");" in the streamer_poisson folder
mesh.f = mesh98.f;
mesh.boundaryexpr = {@(p) abs(p(2,:))<1e-8, @(p) abs(p(1,:)-125)<1e-8, @(p) abs(p(2,:)-125)<1e-8, @(p) abs(p(1,:))<1e-8};
mesh.boundarycondition = [1;2;3;4]; % Set boundary condition for each boundary
mesh.porder = pde.porder;
mesh.dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);


% Initial condition
initu_func_set = {@initu_func_electrons;@initu_func_ions;0};
load 'poissonIC98k.mat';
% load 'poissonIC350k.mat';
UDG_poisson = sol;      % Load in the poisson as UDG
UDG0 = initu_streamer(mesh,initu_func_set,pde.physicsparam);      % Change to UDG0 and UH0 for digaso
UDG0(:,[1,2],:) = UDG0(:,[1,2],:)/ne_star;
UDG0(:,3,:) = UDG_poisson(:,1,:);   
mesh.udg = UDG0;
% UH0=inituhat(master,mesh.elcon,UDG0,app.ncu);
% [QDG, qq, MiCE] = getq(master, mesh, UDG0, UH0, [], 1);
% UDG0(:,app.ncu+1:app.nc,:) = QDG;

% search compilers and set options
pde = setcompilers(pde);       
    
% generate input files and store them in datain folder
[pde,mesh,master,dmd] = preprocessing(pde,mesh);

[sol,pde,mesh] = exasim(pde,mesh);
% return;
% scaplot(mesh,sol(:,1,:),[],2); axis on; axis equal; axis tight;