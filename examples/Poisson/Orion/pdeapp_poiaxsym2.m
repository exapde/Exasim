% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel_poiaxsym2";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
pde.platform = "cpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 4;             % number of MPI processors
pde.hybrid = 1;               % 0 -> LDG, 1 -> HDG
pde.debugmode = 0;
pde.nd = 2;

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;             % polynomial degree
pde.pgauss = 2*pde.porder;
pde.physicsparam = [1 0 1 1];       % unit thermal conductivity
pde.tau = 1.0;              % DG stabilization parameter
pde.linearsolvertol = 1e-8; % GMRES tolerance
pde.ppdegree = 1;          % degree of polynomial preconditioner
pde.RBdim = 0;

[p, t, mesh.dgnodes, XA, YA, XB, YB] = mshOrion(pde.porder, 20, 25, 10, 30);
mesh.p = p';
mesh.t = t';

%mesh = mesh2;
mesh.p(2,:) = mesh.p(2,:) + 1e-4;
mesh.dgnodes(:,2,:) = mesh.dgnodes(:,2,:) + 1e-4;
% expressions for domain boundaries
% x2 = 1e-4, x1 = XB, x1 < 0 | x2 > 3, |x1| < XB
% axis symmetric, outflow, inflow, wall
mesh.boundaryexpr = {@(p) abs(p(2,:)-1e-4)<1e-6, @(p) abs(p(1,:)-XB)< 1e-6, @(p) ((p(1,:) < -1e-6) | (p(2,:) >3)), @(p) abs(p(1,:))< XB + 1e-6};
%mesh.boundaryexpr = {@(p) abs(p(2,:)-1e-4)<1e-6, @(p) abs(p(1,:)-7)< 1e-2, @(p) ((p(1,:) < -1e-6) | (p(2,:) >3)), @(p) abs(p(1,:))< 100 + 1e-6};
mesh.boundarycondition = [1, 4, 3, 2]; % Set boundary condition for each boundary

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh,master,dmd,comstr,runstr] = exasim(pde,mesh);

% plot solution
mesh.porder = pde.porder;
figure(3); clf; scaplot(mesh,sol(:,1,:),[],2); axis on; axis equal; axis tight;

return;

