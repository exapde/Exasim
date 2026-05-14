% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";
pde.modelfile = "pdemodel";

% Choose computing platform and set number of processors
pde.platform = "cpu";
pde.mpiprocs = 1;
pde.hybrid = 1;
pde.debugmode = 0;

E = 1.0;
nu = 0.30;
mu = E / (2 * (1 + nu));
lambda = nu * E / ((1 + nu) * (1 - 2 * nu));
amp = 0.10;

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 2;
pde.pgauss = 2 * pde.porder;
pde.physicsparam = [mu lambda amp];
pde.tau = 2 * mu;
pde.linearsolvertol = 1e-8;
pde.ppdegree = 4;
pde.RBdim = 0;
pde.GMRESrestart = 100;
pde.linearsolveriter = 200;
pde.preconditioner = 1;

mesh = mkmesh_platefluid(pde.porder);
mesh.boundarycondition = [1; 1; 1; 1; 2];

figure(1); clf;
meshplot(mesh, 1);
axis equal;
axis on;
xlabel('x');
ylabel('y');
title('Initial fluid mesh');

[sol,pde,mesh,master,dmd] = exasim(pde,mesh);

mesh.elemtype = 0;
mesh.dgnodes = mesh.dgnodes + sol(:,1:2,:);

figure(2); clf;
meshplot(mesh, 1);
axis equal;
axis on;
xlabel('x');
ylabel('y');
title('Deformed fluid mesh');
