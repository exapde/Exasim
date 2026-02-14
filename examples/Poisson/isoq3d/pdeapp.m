% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";       % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel"; % name of a file defining the PDE model

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 2;             % polynomial degree
% unit thermal conductivity, zero Dirichlet data, and Neumman data
pde.physicsparam = [1 0 1];
pde.tau = 2.0;              % DG stabilization parameter
pde.linearsolveriter = 151;
pde.GMRESrestart = 150;
pde.NLtol = 1e-4;
pde.linearsolvertol = 1e-6;
pde.GMRESortho = 1;
pde.RBdim = 0;
pde.NLiter = 2;
pde.ppdegree=20;

% Choose computing platform and set number of processors
%pde.platform = "gpu";      % choose this option if you want to run the C++ code on Nvidia GPUs
pde.mpiprocs = 4;           % number of MPI processors
pde.hybrid = 1;             % 0 -> LDG, 1-> HDG

mesh = mkmesh_isoq3d(pde.porder, 12);
mesh.boundarycondition = [3;3;3;1;2]; % Set boundary condition for each boundary

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh,master,dmd] = exasim(pde,mesh);

% visualize the numerical solution of the PDE model using Paraview
pde.visscalars = {"temperature", 1};                % list of scalar fields for visualization
pde.visvectors = {"temperature gradient", [2 3 4]}; % list of vector fields for visualization
vis(sol,pde,mesh);                        % visualize the numerical solution
disp("Done!");

