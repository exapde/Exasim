function mesh1 = radaptivity(mesh, rho, drhodx, drhody, diffcoeff)

% initialize pde structure and mesh structure
pde = initializeexasim();

pde.nd = size(mesh.dgnodes, 2);
pde.elemtype = mesh.elemtype;
pde.porder = mesh.porder;      
pde.nodetype = 1;
pde.pgauss = 2*pde.porder;

master = Master(pde);
L = averagevector(master,mesh);
theta = sum(L(:).*rho(:))/sum(L(:));

mesh.xpe = master.xpe;
mesh.telem = master.telem;
mesh.vdg = rho/theta;
mesh.vdg(:,2,:) = drhodx/theta;
mesh.vdg(:,3,:) = drhody/theta;

% solve Poisson equation with  homogeneous Neumann boundary condition
% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel_poisson";    % name of a file defining the PDE model
pde.platform = "cpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;             % number of MPI processors
pde.hybrid = 1;               % 0 -> LDG, 1 -> HDG
pde.debugmode = 0;

% Set discretization parameters, physical parameters, and solver parameters
pde.physicsparam = [1 theta];  
pde.tau = 1.0;              % DG stabilization parameter
pde.linearsolvertol = 1e-8; % GMRES tolerance
pde.preconditioner = 1;
pde.GMRESrestart = 50;

% call exasim to generate and run C++ code to solve the PDE model
[sol, pde] = exasim(pde,mesh);

% get the velocity field as the gradient of the solution of the poisson equation
mesh.vdg(:,4:5,:) = sol(:,2:3,:); 

% mesh1 is the adaptive mesh
mesh1 = mesh; 

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pdet = pde;
pdet.model = "ModelD";          % ModelC, ModelD, ModelW
pdet.modelfile = "pdemodel_transport";    % name of a file defining the PDE model

% Set discretization parameters, physical parameters, and solver parameters
pdet.torder = 2;          % time-stepping order of accuracy
pdet.nstage = 2;          % time-stepping number of stages
pdet.tau = 2.0;               % DG stabilization parameter
pdet.dt = 0.02*ones(1,50);   % time step sizes
pdet.soltime = 1:length(pdet.dt); % steps at which solution are collected
pdet.physicsparam = [diffcoeff theta];  

% solve the transport equation with initial condition u(x,y,t=0) = x
mesh.udg = mesh.dgnodes(:,1,:);
[solt,pdet] = exasim(pdet,mesh);

% the x-component of the adaptive mesh is the solution of the transport
% equation at time t = 1
mesh1.dgnodes(:,1,:) = solt(:,1,:,end);

% solve the transport equation with initial condition u(x,y,t=0) = y
mesh.udg = mesh.dgnodes(:,2,:);
solt = exasim(pdet,mesh);

% the y-component of the adaptive mesh is the solution of the transport
% equation at time t = 1
mesh1.dgnodes(:,2,:) = solt(:,1,:,end);

% plot the velocity field
figure(1); clf; scaplot(mesh,sol(:,2,:),[],2); axis on; axis equal; axis tight;
figure(2); clf; scaplot(mesh,sol(:,3,:),[],2); axis on; axis equal; axis tight;

% plot the mesh density function
figure(3); clf; scaplot(mesh,rho(:,1,:),[],1); axis on; axis equal; axis tight;

% plot the original mesh
figure(4); clf; meshplot(mesh,1); axis on; axis equal; axis tight;

% plot the adaptive mesh
figure(5); clf; meshplot(mesh1,1); axis on; axis equal; axis tight;

