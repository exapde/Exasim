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
pde.porder = 3;             % polynomial degree
pde.pgauss = 2*pde.porder;  % gauss quad order
pde.physicsparam = [1 1 0 0 0.01]; 
pde.tau = 2.0;              % DG stabilization parameter
pde.RBdim = 0;              % reduced basis dimension for preconditioner
pde.linearsolvertol = 1e-6; % GMRES tolerance
pde.NLtol = 1e-8; % GMRES tolerance
pde.GMRESortho = 1;
pde.GMRESrestart = 100;
pde.linearsolveriter = 1000;
pde.preconditioner = 1;
pde.ppdegree = 0;           % degree of polynomial preconditioner

% create a grid of 8 by 8 by 8 hexes on the unit cube
[mesh.p,mesh.t] = cubemesh(20,4,4,1);
mesh.p(1,:) = 5*mesh.p(1,:);

% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(2,:))<1e-6, @(p) abs(p(1,:)-5)<1e-6, @(p) abs(p(2,:)-1)<1e-6, @(p) abs(p(1,:))<1e-6, @(p) abs(p(3,:))<1e-6, @(p) abs(p(3,:)-1)<1e-6};
mesh.boundarycondition = [1;1;1;2;1;3]; % Set boundary condition for each boundary

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh,master,dmd] = exasim(pde,mesh);

bndexpr = {'all(p(:,1)<1e-6)','all(p(:,1)>max(p0(:,1))-1e-6)', ...
           'all(p(:,2)<1e-6)','all(p(:,2)>max(p0(:,2))-1e-6)', ...
           'all(p(:,3)<1e-6)','all(p(:,3)>max(p0(:,3))-1e-6)'};     
mesh1 = mkmesh(mesh.p',mesh.t',pde.porder,bndexpr,1,1);
mesh1.p = mesh1.p';
mesh1.t = mesh1.t';
figure(1);clf;meshplot(mesh1,1); axis on;
mesh1.dgnodes = mesh1.dgnodes + sol(:,1:3,:); 
figure(2);clf;meshplot(mesh1,1); 
axis equal; axis on; box on; 
set(gca,'fontsize', 16);

xlabel('x');
ylabel('y');
zlabel('z');

