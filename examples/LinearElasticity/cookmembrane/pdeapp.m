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
pde.mpiprocs = 1;           % number of MPI processors
pde.hybrid = 1;             % 0 -> LDG, 1-> HDG
pde.debugmode = 0;

E=240.565;
nu=0.495;
mu=E/(2*(1+nu));
lambda=nu*E/((1+nu)*(1-2*nu));

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;             % polynomial degree
pde.pgauss = 2*pde.porder;  % gauss quad order
pde.physicsparam = [mu lambda 0 16]; 
pde.tau = 2*mu;              % DG stabilization parameter
pde.linearsolvertol = 1e-6; % GMRES tolerance
pde.ppdegree = 10;           % degree of polynomial preconditioner
pde.RBdim = 0;              % reduced basis dimension for preconditioner
pde.GMRESrestart = 100;
pde.linearsolveriter = 100;
pde.preconditioner = 1;

% create a grid of 8 by 8 by 8 hexes on the unit cube
elemtype = 1;
[mesh.p,mesh.t] = cookmembranegrid(10,10,elemtype);
xmin = min(mesh.p(1,:));
xmax = max(mesh.p(1,:));

% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(1,:)-xmin)<1e-6, @(p) abs(p(1,:)-xmax)<1e-6, @(p) p(2,:)<44+1e-6, @(p) p(2,:)>44-1e-6};
mesh.boundarycondition = [2;3;1;1]; % Set boundary condition for each boundary

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh,master,dmd] = exasim(pde,mesh);

mesh.elemtype = elemtype;
figure(1);clf;meshplot(mesh,1); axis on;
mesh.dgnodes = mesh.dgnodes + sol(:,1:2,:); 
figure(2);clf;meshplot(mesh,1); axis on;
xlabel('x');
ylabel('y');

