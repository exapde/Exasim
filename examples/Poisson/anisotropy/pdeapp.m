% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
pde.platform = "cpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;             % number of MPI processors
pde.hybrid = 1;               % 0 -> LDG, 1 -> HDG
pde.debugmode = 0;
pde.nd = 2;

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;             % polynomial degree
pde.pgauss = 2*pde.porder;
pde.preconditioner = 1;
pde.ppdegree = 1;          % degree of polynomial preconditioner
pde.RBdim = 0;
pde.GMRESrestart = 50;
pde.saveSolBouFreq = 1;
pde.ibs = 1;

% create a grid of 8 by 8 on the unit square
[mesh.p,mesh.t] = squaremesh(20,20,1,0);

L = 10;
pde.physicsparam = [1, 1];       % unit thermal conductivity
pde.linearsolvertol = 1e-7/L; % GMRES tolerance
mesh.p(1,:) = L*mesh.p(1,:);

pde.tau = 1;              % DG stabilization parameter

% mesht = mkmesh_square(8,8,pde.porder,0,1,1,1,1);
% mesh.dgnodes = mesht.dgnodes;
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(2,:))<1e-8, @(p) abs(p(1,:)-L)<1e-8, @(p) abs(p(2,:)-1)<1e-8, @(p) abs(p(1,:))<1e-8};
mesh.boundarycondition = [1;1;1;1]; % Set boundary condition for each boundary

pde.gencode = 0;
% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh,master,dmd] = exasim(pde,mesh);

% plot solution
mesh.porder = pde.porder;
mesh.dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);
figure(1); clf; scaplot(mesh,sol(:,1,:),[],2); axis on; axis equal; axis tight;

mesh0 = mesh;
sol0 = sol;

