% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pdehm structure and mesh structure
[pdehm,~] = initializeexasim();
pdehm.buildpath = string(pwd()) + "/hm";

% Define a pdehm model: governing equations, initial solutions, and boundary conditions
pdehm.model = "ModelD";          % ModelC, ModelD, ModelW
pdehm.modelfile = "pdemodel_hm";    % name of a file defining the pdehm model

% Choose computing platform and set nuomber of processors
pdehm.platform = "cpu";         % choose this option if NVIDIA GPUs are available
pdehm.mpiprocs = 1;             % number of MPI processors
pdehm.hybrid = 1;               % 0 -> LDG, 1 -> HDG
pdehm.debugmode = 0;
pdehm.nd = 2;

% Set discretization parameters, physical parameters, and solver parameters
pdehm.porder = 4;             % polynomial degree
pdehm.pgauss = 2*pdehm.porder;
pdehm.physicsparam = kappa(1)^2*3e-3;       
pdehm.tau = 1.0;              % DG stabilization parameter
pdehm.linearsolvertol = 1e-8; % GMRES tolerance
pdehm.ppdegree = 1;          % degree of polynomial preconditioner
pdehm.RBdim = 0;

meshhm = mkmesh_cyl(pdehm.porder);
meshhm.boundarycondition = [1;1;1]; 
div = divergence(sol, 1);
meshhm.vdg = limiting(div,0,3,1e3,0);

% call exasim to generate and run C++ code to solve the pdehm model
[solhm,pdehm,meshhm] = exasim(pdehm,meshhm);
s = solhm(:,1,:);
s = s/max(s(:));
av = (s-S0).*(atan(gamma*(s-S0))/pi + 0.5) - atan(gamma)/pi + 0.5;    
dist = tanh(mesh.dist*5);
av = lambda(1)*(av.*dist);         

% plot solution
figure(2); clf; scaplot(meshhm,av); axis on; axis equal; axis tight;

