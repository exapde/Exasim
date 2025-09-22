% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pdeht structure and mesh structure
[pdeht,~] = initializeexasim();
pdeht.buildpath = string(pwd()) + "/ht";

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pdeht.model = "ModelD";          % ModelC, ModelD, ModelW
pdeht.modelfile = "pdemodel_axialht";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
pdeht.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pdeht.mpiprocs = 1;             % number of MPI processors
pdeht.hybrid = 1;               % 0 -> LDG, 1 -> HDG
pdeht.debugmode = 0;
pdeht.nd = 2;
pdeht.gencode = 1;

kappa = 0.2;

% Set discretization parameters, physical parameters, and solver parameters
pdeht.porder = 3;             % polynomial degree
pdeht.pgauss = 2*pdeht.porder;
pdeht.physicsparam = [kappa Twall/Tref * Tinf];       % unit thermal conductivity
pdeht.tau = 1.0;              % DG stabilization parameter
pde.GMRESrestart = 100;     
pde.linearsolvertol = 1e-8;   % GMRES tolerance
pde.linearsolveriter = 100; 
pdeht.ppdegree = 10;          % degree of polynomial preconditioner
pdeht.RBdim = 0;

[~, meshht] = mkmesh_orion(pdeht.porder);
meshht.p(2,:) = meshht.p(2,:) + 1e-3;
meshht.dgnodes(:,2,:) = meshht.dgnodes(:,2,:) + 1e-3;
meshht.boundaryexpr = {@(p) abs(p(2,:)-1e-3)<1e-6, @outerwall, @(p) abs(p(1,:))< 20 + 1e-6};
% axis symmetric (homogenuous Neumann), heat flux (nonhomogenuous Neumann), thermal (nonhomogenuous Dirichlet) 
meshht.boundarycondition = [3, 4, 1]; % Set boundary condition for each boundary

return;





