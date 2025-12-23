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

% Set discretization parameters, physical parameters, and solver parameters
pde.dae_alpha = 0.0;
pde.dae_beta = 0.0;
pde.porder = 3;             % polynomial degree
pde.pgauss = 2*pde.porder;  % gauss quad order
pde.tau = 100.0;              % DG stabilization parameter
pde.linearsolvertol = 1e-6; % GMRES tolerance
pde.ppdegree = 10;           % degree of polynomial preconditioner
pde.RBdim = 0;              % reduced basis dimension for preconditioner
pde.GMRESrestart = 100;
pde.linearsolveriter = 100;
pde.preconditioner = 1;
pde.nd = 2;

% create a grid of 8 by 8 by 8 hexes on the unit cube
pde.elemtype = 1;
[mesh.p,mesh.t] = squaremesh(50,50,1,pde.elemtype);
mesh.p(1,:) = 2*mesh.p(1,:) - 1;
mesh.p(2,:) = 2*mesh.p(2,:) - 1;

% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(2,:)+1)<1e-8, @(p) abs(p(1,:)-1)<1e-8, @(p) abs(p(2,:)-1)<1e-8, @(p) abs(p(1,:)+1)<1e-8};
mesh.boundarycondition = [1;1;1;1]; % Set boundary condition for each boundary
mesh.porder = pde.porder;
mesh.f = facenumbering(mesh.p,mesh.t,pde.elemtype,mesh.boundaryexpr,mesh.periodicexpr);
mesh.dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);

master = Master(pde);

x = (mesh.dgnodes(:,1,:));
y = (mesh.dgnodes(:,2,:));
r = sqrt(x.^2+y.^2);
a1 = 20;
a2 = 100;
a = 0.25;
rho = 1 + a1*sech(a2*(r.^2-a^2));
L = averagevector(master,mesh);
theta = sum(L(:).*rho(:))/sum(L(:));
pde.physicsparam = [0.025 1 theta a1 a2 a];  

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh,master,dmd] = exasim(pde,mesh);

mesh.elemtype = pde.elemtype;
figure(6); clf; scaplot(mesh, sol(:,1,:),[],2); xlabel('x'); ylabel('y');
figure(7); clf; scaplot(mesh, sol(:,2,:),[],2); xlabel('x'); ylabel('y');


