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
pde.physicsparam = [0.1 1e4]; 
pde.tau = 2.0;              % DG stabilization parameter
pde.linearsolvertol = 1e-6; % GMRES tolerance
pde.ppdegree = 10;           % degree of polynomial preconditioner
pde.RBdim = 0;              % reduced basis dimension for preconditioner
pde.GMRESrestart = 100;
pde.linearsolveriter = 100;
pde.preconditioner = 1;

% create a grid of 8 by 8 by 8 hexes on the unit cube
elemtype = 1;
[mesh.p,mesh.t] = squaremesh(16,16,1,elemtype);

% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(2,:))<1e-6, @(p) abs(p(1,:)-1)<1e-6, @(p) abs(p(2,:)-1)<1e-6, @(p) abs(p(1,:))<1e-6};
mesh.boundarycondition = [1;1;1;1]; % Set boundary condition for each boundary

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh,master,dmd] = exasim(pde,mesh);

mesh.elemtype = elemtype;
figure(1); clf; scaplot(mesh, sol(:,1,:),[],2); xlabel('x'); ylabel('y');
figure(2); clf; scaplot(mesh, sol(:,2,:),[],2); xlabel('x'); ylabel('y');

x1 = mesh.dgnodes(:,1,:);
x2 = mesh.dgnodes(:,2,:);
Re = 1/pde.physicsparam(1); 
lam = Re/2 - sqrt(Re^2/4+4*pi^2);
u1 = 1-exp(lam*x1).*cos(2*pi*x2);
u2 = (lam/(2*pi))*exp(lam*x1).*sin(2*pi*x2);
figure(3); clf; scaplot(mesh,u1,[],2); xlabel('x'); ylabel('y');
figure(4); clf; scaplot(mesh,u2,[],2); xlabel('x'); ylabel('y');

