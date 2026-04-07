% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";       % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel"; % name of a file defining the PDE model

% Choose computing platform and set number of processors
pde.platform = "cpu";       % choose this option if you want to run the C++ code on Nvidia GPUs
pde.mpiprocs = 8;           % number of MPI processors
pde.hybrid = 1;             % 0 -> LDG, 1-> HDG
pde.debugmode = 0;

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 2;             % polynomial degree
pde.pgauss = 2*pde.porder;  % gauss quad order
pde.physicsparam = [1.25 0.833 0.015 0 -0.015]; 
pde.tau = 5.0;              % DG stabilization parameter
pde.RBdim = 0;              % reduced basis dimension for preconditioner
pde.linearsolvertol = 1e-6; % GMRES tolerance
pde.NLtol = 1e-6; % GMRES tolerance
pde.GMRESortho = 1;
pde.GMRESrestart = 100;
pde.linearsolveriter = 200;
pde.preconditioner = 1;
pde.ppdegree = 0;           % degree of polynomial preconditioner

% create a grid of 8 by 8 by 8 hexes on the unit cube
mesh = mkmesh_seawall3d3(pde.porder);
mesh.boundarycondition = [2;2;2;1;3;1;1;2]; % Set boundary condition for each boundary

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh,master,dmd] = exasim(pde,mesh);

pde.visscalars = {"x displacement", 1};                % list of scalar fields for visualization
pde.visvectors = {"displacement field", [1 2 3]}; % list of vector fields for visualization
mesh.xpe = master.xpe;
mesh1 = mesh;
mesh1.dgnodes = mesh.dgnodes + sol(:,1:3,:);
dgnodes = vis(sol,pde,mesh1);  % visualize the numerical solution

