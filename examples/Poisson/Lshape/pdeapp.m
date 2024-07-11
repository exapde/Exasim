% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Install/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();  

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
%pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 4;              % number of MPI processors
pde.hybrid = 1;                % 0 -> LDG, 1-> HDG
pde.ppdegree = 10;             % degree of polynomial preconditioner

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;          % polynomial degree
pde.physicsparam = 1;    % unit thermal conductivity
pde.tau = 1.0;           % DG stabilization parameter

% call Gmsh to generate a mesh on L-shaped domain, see lshape.geo for details
dim = 2; elemtype = 0;
[mesh.p, mesh.t] = gmshcall(pde, "lshape", dim, elemtype);
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(1,:))<2};
mesh.boundarycondition = [1]; % Set boundary condition for each boundary

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh] = exasim(pde,mesh);

% visualize the numerical solution of the PDE model using Paraview
% pde.visscalars = {"temperature", 1};  % list of scalar fields for visualization
% pde.visvectors = {"temperature gradient", [2 3]}; % list of vector fields for visualization
% vis(sol,pde,mesh); % visualize the numerical solution
% disp("Done!");

% plot solution
mesh.porder = pde.porder;
mesh.dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);
figure(1); clf; scaplot(mesh,sol(:,1,:),[],2,1); axis on; axis equal; axis tight;


        
