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
pdeht.platform = "cpu";         % choose this option if NVIDIA GPUs are available
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
pdeht.ppdehtgree = 1;          % degree of polynomial preconditioner
pdeht.RBdim = 0;

meshht = mkmesh_cylht(pdeht.porder);
meshht.boundarycondition = [4;3;1];
meshht.vdg = zeros(size(meshht.dgnodes,1),4,size(meshht.dgnodes,3));
meshht.vdg(:,1,:) = 0.01;

% call exasim to generate and run C++ code to solve the PDE model
[solht,pdeht,meshht,masterht,dmdht] = exasim(pdeht,meshht);

% plot solution
figure(2); clf; scaplot(meshht, solht(:,1,:));

return;





