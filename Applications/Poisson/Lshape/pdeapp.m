% Specify an Exasim version to run
version = "Version0.1";

% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim(version);  

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;          % polynomial degree
pde.physicsparam = 1;    % unit thermal conductivity
pde.tau = 1.0;           % DG stabilization parameter

% Choose computing platform and set number of processors
%pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors

% call Gmsh to generate a mesh on L-shaped domain, see lshape.geo for details
dim = 2; elemtype = 0;
[mesh.p, mesh.t] = gmshcall(pde, "lshape", dim, elemtype);
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(1,:))<2};
mesh.boundarycondition = [1]; % Set boundary condition for each boundary

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh] = exasim(pde,mesh);

% visualize the numerical solution of the PDE model using Paraview
pde.visscalars = {"temperature", 1};  % list of scalar fields for visualization
pde.visvectors = {"temperature gradient", [2 3]}; % list of vector fields for visualization
vis(sol,pde,mesh); % visualize the numerical solution
disp("Done!");

        
