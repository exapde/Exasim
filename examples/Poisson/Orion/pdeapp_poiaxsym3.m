% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel_poiaxsym1";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
pde.platform = "cpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;             % number of MPI processors
pde.hybrid = 1;               % 0 -> LDG, 1 -> HDG
pde.debugmode = 0;
pde.nd = 2;

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;             % polynomial degree
pde.pgauss = 2*pde.porder;
pde.physicsparam = [1 0 1 1];       % unit thermal conductivity
pde.tau = 1.0;              % DG stabilization parameter
pde.linearsolvertol = 1e-8; % GMRES tolerance
pde.ppdegree = 10;          % degree of polynomial preconditioner
pde.GMRESrestart = 50;
pde.RBdim = 0;

% % create a grid of 8 by 8 on the unit square
% [p, t, mesh.dgnodes, XA, YA, XB, YB] = mshOrion(pde.porder, 20, 25, 10, 30);
% mesh.p = p';
% mesh.t = t';
% % expressions for domain boundaries
% % x2 = 0, x1 = XB, x1 < 0 | x2 > 3, |x1| < XB
% % axis symmetric, outflow, inflow, wall
% mesh.boundaryexpr = {@(p) abs(p(2,:))<1e-6, @(p) abs(p(1,:)-XB)< 1e-6, @(p) ((p(1,:) < 1e-6) | (p(2,:) >3)), @(p) abs(p(1,:))< XB + 1e-6};
% mesh.boundarycondition = [1, 4, 3, 2]; % Set boundary condition for each boundary

% [~, mesh] = mkmesh_orion5(pde.porder, 0.12);
% mesh.boundarycondition = [1, 3, 2];

%[mesh, ~] = mkmesh_orion(pde.porder, 0.12);
mesh = mkmesh_orion6(pde.porder);
mesh.boundarycondition = [1, 4, 3, 2];

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh,master,dmd,comstr,runstr] = exasim(pde,mesh);

% plot solution
mesh.porder = pde.porder;
figure(2); clf; scaplot(mesh,sol(:,1,:),[],2); axis on; axis equal; axis tight;

return;

