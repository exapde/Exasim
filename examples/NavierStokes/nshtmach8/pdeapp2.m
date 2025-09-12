% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde{2} structure and mesh structure
[pde{2},~] = initializeexasim();
pde{2}.buildpath = string(pwd());

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde{2}.model = "ModelD";          % ModelC, ModelD, ModelW
pde{2}.modelfile = "pdemodel2";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
pde{2}.platform = "cpu";         % choose this option if NVIDIA GPUs are available
pde{2}.mpiprocs = 1;             % number of MPI processors
pde{2}.hybrid = 1;               % 0 -> LDG, 1 -> HDG
pde{2}.debugmode = 0;
pde{2}.nd = 2;
pde{2}.gencode = 1;

kappa = 0.2;

% Set discretization parameters, physical parameters, and solver parameters
pde{2}.porder = 4;             % polynomial degree
pde{2}.pgauss = 2*pde{2}.porder;
pde{2}.physicsparam = [kappa Twall/Tref * Tinf];       % unit thermal conductivity
pde{2}.tau = 1.0;              % DG stabilization parameter
pde{2}.GMRESrestart = 100;         %try 50
pde{2}.linearsolvertol = 1e-8; % GMRES tolerance
pde{2}.linearsolveriter = 100; %try 100
pde{2}.RBdim = 0;
pde{2}.ppdegree = 20;
pde{2}.NLtol = 1e-6;              % Newton tolerance
pde{2}.NLiter = 30;                 % Newton iterations
pde{2}.matvectol=1e-6;             % tolerance for matrix-vector multiplication

mesh{2} = mkmesh_cylht(pde{2}.porder);
mesh{2}.udg = zeros(size(mesh{2}.dgnodes,1),3,size(mesh{2}.dgnodes,3));
mesh{2}.udg(:,1,:) = Twall/Tref * Tinf;





