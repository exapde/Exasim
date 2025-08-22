% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel_poisson";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
pde.platform = "cpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;             % number of MPI processors
pde.hybrid = 1;               % 0 -> LDG, 1 -> HDG
pde.debugmode = 0;

% Set discretization parameters, physical parameters, and solver parameters
pde.tau = 1.0;              % DG stabilization parameter
pde.linearsolvertol = 1e-8; % GMRES tolerance
pde.preconditioner = 1;
pde.GMRESrestart = 50;
pde.ppdegree = 1;          % degree of polynomial preconditioner
pde.RBdim = 0;

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh,master,dmd,comstr,runstr] = exasim(pde,mesh);
