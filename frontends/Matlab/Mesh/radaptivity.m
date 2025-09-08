function mesh1 = radaptivity(mesh, rho, drhodx, drhody, diffcoeff)

if nargin<5
    diffcoeff = 0;
end

% initialize pde structure and mesh structure
[pde, mesh0] = initializeexasim();
rpath = pde.buildpath + "/radaptivity";
if exist(rpath,'dir') == 0
    mkdir(rpath);
end

% get only necessary fields
mesh0.p = mesh.p;
mesh0.t = mesh.t;
mesh0.boundaryexpr = mesh.boundaryexpr; 
mesh0.boundarycondition = ones(length(mesh.boundarycondition),1); 
mesh0.porder = mesh.porder;
mesh0.elemtype = mesh.elemtype;
mesh0.f = mesh.f;
mesh0.dgnodes = mesh.dgnodes;
mesh = mesh0;

pde.nd = size(mesh.dgnodes, 2);
pde.elemtype = mesh.elemtype;
pde.porder = mesh.porder;      
pde.nodetype = 1;
pde.pgauss = 2*pde.porder;

pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel_poisson";    % name of a file defining the PDE model
pde.platform = "cpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;             % number of MPI processors
pde.hybrid = 1;               % 0 -> LDG, 1 -> HDG
pde.debugmode = 0;

master = Master(pde);
mesh.xpe = master.xpe;
mesh.telem = master.telem;

L = averagevector(master,mesh);
theta = sum(L(:).*rho(:))/sum(L(:));
mesh.vdg = rho/theta;
mesh.vdg(:,2,:) = drhodx/theta;
mesh.vdg(:,3,:) = drhody/theta;

% Set discretization parameters, physical parameters, and solver parameters
pde.physicsparam = [1 theta];  
pde.tau = 2.0;              % DG stabilization parameter
pde.linearsolvertol = 1e-8; % GMRES tolerance
pde.preconditioner = 1;
pde.GMRESrestart = 200;
pde.linearsolveriter = 400;
pde.ppdegree = 10;
pde.NLtol = 1e-8;

% solve Poisson equation with  homogeneous Neumann boundary condition
pde.buildpath = rpath + "/poisson"; 
[pde,mesh,master,dmd] = preprocessing(pde,mesh);
if exist(pde.buildpath + "/cpuEXASIM", "file") == 0
  kkgencode(pde);
  cmakecompile(pde); 
end           
runcode(pde, 1); % run C++ code
% get solution from output files in dataout folder
sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');

% plot the velocity field
figure(1); clf; scaplot(mesh,sol(:,2,:),[],2); axis on; axis equal; 
figure(2); clf; scaplot(mesh,sol(:,3,:),[],2); axis on; axis equal; 

[~,cgelcon,rowent2elem,colent2elem] = mkcgent2dgent(mesh.dgnodes,1e-8);
[~, ~, jac] = volgeom(master.shapent,permute(mesh.dgnodes,[1 3 2]));
jac = reshape(jac,[],1,size(mesh.dgnodes,3));
jac = dg2cg2(jac, cgelcon, colent2elem, rowent2elem);
jac = dg2cg2(jac, cgelcon, colent2elem, rowent2elem);

% get the velocity field as the gradient of the solution of the poisson equation
mesh.vdg(:,4:5,:) = sol(:,2:3,:); 
mesh.vdg(:,6,:) = reshape(jac, size(mesh.dgnodes(:,1,:)));

% mesh1 is the adaptive mesh
mesh1 = mesh; 

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pdet = pde;
pdet.model = "ModelD";          % ModelC, ModelD, ModelW
pdet.modelfile = "pdemodel_transport";    % name of a file defining the PDE model

% Set discretization parameters, physical parameters, and solver parameters
pdet.torder = 2;          % time-stepping order of accuracy
pdet.nstage = 2;          % time-stepping number of stages
pdet.tau = 1.0;               % DG stabilization parameter
pdet.dt = 0.02*ones(50,1);   % time step sizes
pdet.soltime = 1:length(pdet.dt); % steps at which solution are collected
pdet.physicsparam = [diffcoeff theta];  
pdet.RBdim = 1;
pdet.ppdegree = 1;

% solve the transport equation with initial condition u(x,y,t=0) = x
%mesh.boundarycondition(3:5) = 2;
mesh.udg = mesh.dgnodes(:,1,:);
%[solt,pdet] = exasim(pdet,mesh);
pdet.buildpath = rpath + "/transport"; 
[pdet,mesh,master,dmd] = preprocessing(pdet,mesh);
if exist(pdet.buildpath + "/cpuEXASIM", "file") == 0
  kkgencode(pdet);
  cmakecompile(pdet); 
end           
runcode(pdet, 1); % run C++ code
% get solution from output files in dataout folder
solt = fetchsolution(pdet,master,dmd, pdet.buildpath + '/dataout');

% the x-component of the adaptive mesh is the solution of the transport
% equation at time t = 1
mesh1.dgnodes(:,1,:) = solt(:,1,:,end);

% solve the transport equation with initial condition u(x,y,t=0) = y
%mesh.boundarycondition(3:5) = 3;
mesh.udg = mesh.dgnodes(:,2,:);
%solt = exasim(pdet,mesh);
[pdet,mesh,master,dmd] = preprocessing(pdet,mesh);
runcode(pdet, 1); % run C++ code
solt = fetchsolution(pdet,master,dmd, pdet.buildpath + '/dataout');

% the y-component of the adaptive mesh is the solution of the transport
% equation at time t = 1
mesh1.dgnodes(:,2,:) = solt(:,1,:,end);

% plot the mesh density function
figure(3); clf; scaplot(mesh,rho(:,1,:),[],1); axis on; axis equal; 

% plot the original mesh
figure(4); clf; meshplot(mesh,1); axis on; axis equal; 

% plot the adaptive mesh
figure(5); clf; meshplot(mesh1,1); axis on; axis equal; 

