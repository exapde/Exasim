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
pde.physicsparam = [1 0 1 1];       % unit thermal conductivity
pde.tau = 1.0;              % DG stabilization parameter
pde.linearsolvertol = 1e-8; % GMRES tolerance
pde.ppdegree = 10;          % degree of polynomial preconditioner
pde.GMRESrestart = 50;
pde.RBdim = 0;

mesh = mkmesh_sharpb2(pde.porder);
mesh.boundarycondition = [1, 3, 2, 4];

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh,master,dmd,comstr,runstr] = exasim(pde,mesh);

% plot solution
mesh.porder = pde.porder;
figure(2); clf; scaplot(mesh,sol(:,1,:),[],2); axis on; axis equal; axis tight; colorbar;

return;

