% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";       % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel"; % name of a file defining the PDE model

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 2;             % polynomial degree
% unit thermal conductivity, zero Dirichlet data, and Neumman data
pde.physicsparam = [1 0.0 1.0 0.0 0.0];
pde.tau = 2.0;              % DG stabilization parameter
pde.linearsolveriter = 400;
pde.GMRESrestart = 200;
pde.NLtol = 1e-4;
pde.GMRESortho = 1;
pde.RBdim = 1;
pde.NLiter = 8;

% Choose computing platform and set number of processors
%pde.platform = "gpu";      % choose this option if you want to run the C++ code on Nvidia GPUs
pde.mpiprocs = 4;           % number of MPI processors

% use Gmsh to generate a mesh
nd = 3; elemtype = 0;
[mesh.p,mesh.t] = gmshcall(pde, "coneincube", nd, elemtype);
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(2,:)+10)<1e-3, @(p) abs(p(1,:)-10)<1e-3, @(p) abs(p(2,:)-10)<1e-3, @(p) abs(p(1,:)+10)<1e-3, @(p) abs(p(3,:)+10)<1e-3, @(p) abs(p(3,:)-10)<1e-3, @(p) abs(p(3,:))<1e3};
mesh.boundarycondition = [2;2;2;2;2;2;1]; % Set boundary condition for each boundary

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh,master,dmd,compilerstr,runstr] = exasim(pde,mesh);

