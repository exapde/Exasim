% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pdeht structure and mesh structure
[pdeht,~] = initializeexasim();
pdeht.buildpath = string(pwd()) + "/ht";

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pdeht.model = "ModelD";          % ModelC, ModelD, ModelW
pdeht.modelfile = "pdemodel_ht";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
pdeht.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pdeht.mpiprocs = 1;             % number of MPI processors
pdeht.hybrid = 1;               % 0 -> LDG, 1 -> HDG
pdeht.debugmode = 0;
pdeht.nd = 2;
pdeht.gencode = 1;

kappa = 0.2;

% Set discretization parameters, physical parameters, and solver parameters
pdeht.porder = 4;             % polynomial degree
pdeht.pgauss = 2*pdeht.porder;
pdeht.physicsparam = [kappa Twall/Tref * Tinf];       % unit thermal conductivity
pdeht.tau = 1.0;              % DG stabilization parameter
pdeht.linearsolvertol = 1e-8; % GMRES tolerance
pdeht.ppdegree = 1;          % degree of polynomial preconditioner
pdeht.RBdim = 0;

meshht = mkmesh_cylht2(pdeht.porder);
meshht.bndexpr = {'all(sqrt(p(:,1).^2+p(:,2).^2)<1+1e-4)'    'all(p(:,1)>-1e-3)'    'all(abs(p(:,1))<20)'};

meshht.boundarycondition = [4;3;1];
meshht.vdg = zeros(size(meshht.dgnodes,1),4,size(meshht.dgnodes,3));
meshht.vdg(:,1,:) = 0.01;

%%% Get working on adapted mesh
tmp = load("meshht_adapt.mat", "meshht_adapt");
mesh_mat = tmp.meshht_adapt;
[pdeht,meshht2,masterht,~] = preprocessing(pdeht,meshht);

pht = dgnodes_to_p(mesh_mat.dgnodes(:,1:2,:), meshht, masterht);
tht = mesh_mat.t;

meshht2 = mkmesh(pht,tht,porder,meshht.bndexpr,1,1);
meshht2.p = meshht2.p';
meshht2.t = meshht2.t';
meshht2.dgnodes = mesh_mat.dgnodes(:,1:2,:);

meshht2.porder = porder;
meshht2.boundaryexpr = meshht.boundaryexpr; %TODO
meshht2.bndexpr = meshht.bndexpr;           %TODO
meshht2.periodicexpr = {};

meshht2.boundarycondition = meshht.boundarycondition; 

[pdeht,meshht2,masterht,~] = preprocessing(pdeht,meshht2);
meshht2.vdg = meshht.vdg;

% call exasim to generate and run C++ code to solve the PDE model
[solht,pdeht,meshht,masterht,dmdht] = exasim(pdeht,meshht2);

% plot solution
figure(2); clf; scaplot(meshht, solht(:,1,:));

return;





