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
pde.nd = 1;

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 2;             % polynomial degree
pde.pgauss = 2*pde.porder;
pde.physicsparam = 1;       % unit thermal conductivity
pde.tau = 1.0;              % DG stabilization parameter
pde.linearsolvertol = 1e-8; % GMRES tolerance
pde.ppdegree = 1;          % degree of polynomial preconditioner
pde.RBdim = 0;

% create a 1D grid from 0 to 1
[mesh.p,mesh.t] = linemesh(5);

% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(1))<1e-8, @(p) abs(p(1)-1<1e-8)};

mesh.boundarycondition = [1;1]; %set boundary condition for each boundary

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh,master,dmd,comstr,runstr] = exasim(pde,mesh);
% sol: num pts per element x num components x num elements
u = sol(:,1,:);
figure(1); clf; plot(mesh.dgnodes(:),u(:))
grid on

return;


