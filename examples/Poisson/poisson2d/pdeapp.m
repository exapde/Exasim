% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
pde.platform = "cpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors
pde.hybrid = 1;
pde.debugmode = 0;
pde.nd = 2;

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 2;             % polynomial degree
pde.pgauss = 2*pde.porder;
pde.physicsparam = 1;       % unit thermal conductivity
pde.tau = 1.0;              % DG stabilization parameter
pde.linearsolvertol = 1e-8; % GMRES tolerance
pde.ppdegree = 1;          % degree of polynomial preconditioner

% create a grid of 8 by 8 on the unit square
[mesh.p,mesh.t] = squaremesh(4,4,1,1);
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(2,:))<1e-8, @(p) abs(p(1,:)-1)<1e-8, @(p) abs(p(2,:)-1)<1e-8, @(p) abs(p(1,:))<1e-8};
mesh.boundarycondition = [1;1;1;1]; % Set boundary condition for each boundary

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh] = exasim(pde,mesh);

% plot solution
% mesh.porder = pde.porder;
% figure(1); clf; scaplot(mesh,sol(:,1,:),[],2); axis on; axis equal; axis tight;


% visualize the numerical solution of the PDE model using Paraview
% pde.visscalars = {"temperature", 1};  % list of scalar fields for visualization
% pde.visvectors = {"temperature gradient", [2 3]}; % list of vector fields for visualization
% xdg = vis(sol,pde,mesh); % visualize the numerical solution
% disp("Done!");

master1 = Master(pde);
mesh1 = hdgmesh(mesh, pde.porder);
x = mesh1.dgnodes(:,1,:);
y = mesh1.dgnodes(:,2,:);
u = x.*x + y.*y;
UDG = u;
UDG(:,2,:) = -2*x;
UDG(:,3,:) = -2*y;
UDG = 0*UDG;
UH = inituhat(master1,mesh1.elcon,UDG,1);
pde.bcm = [1; 1; 1; 1];
pde.bcs = [1; 1; 1; 1]*0;
pde.debugmode=1;
[UDG,UH]= hdgsolve(master1,mesh1,pde,UDG,UH,[]);
% compareexasim(master1, mesh1, pde);


% nge = 25;
% [Mass, Minv, C, E, udg] = qEquationExasim(pde, mesh);
% [Ru, B, D, fg, fg_udg, sg, sg_udg] = uEquationElemExasim(pde, mesh, nge);
% 
% run('/Users/ngoccuongnguyen/Dropbox (MIT)/HDGv1.0/setup.m');
% ngrid = 3;
% mesh1   = mkmesh_square(ngrid+1,ngrid+1,pde.porder,0,1,1,1,1);
% master = mkmaster(mesh1,2*pde.porder);
% app1.hybrid = "hdg";
% app1.bcm = [1; 1; 1; 1];
% app1.bcs = [1; 1; 1; 1]*0;
% app1.flg_q = 1;
% app1.flg_p = 1;
% app1.appname = "poisson";
% app1.dt = [];
% app1.nc = 3;
% app1.ncu = 1;
% app1.nch = 1;
% app1.time = 0;
% app1.tdep = 0;
% app1.wave=0;
% app1.alag=0;
% app1.adjoint=0;
% app1.linearproblem=0;
% app1.fc_u = 1;
% app1.fc_q = 1;
% app.denseblock=0;
% app1.arg = {};
% app1.flux = "flux";
% app1.source = "source";
% app1.fbou   = "fbou";
% app1.fhat   = "fhat";
% [master,mesh1] = preprocess(master,mesh1,app1.hybrid);
% [Mass1, Minv1, C1, E1] = qequation(master, mesh1);
% 
% x = mesh1.dgnodes(:,1,:);
% y = mesh1.dgnodes(:,2,:);
% u = x.*x + y.*y;
% UDG = u;
% UDG(:,2,:) = -2*x;
% UDG(:,3,:) = -2*y;
% SDG = [];
% 
% [Ru1, B1, D1] = uequationelem(master,app1,mesh1.dgnodes,UDG,SDG);
% 
% max(abs(Minv(:)-Minv1(:)))
% max(abs(C(:)-C1(:)))
% max(abs(E(:)-E1(:)))
% max(abs(udg(:)-UDG(:)))
% 
% max(abs(Ru(:)-Ru1(:)))
% max(abs(B(:)-B1(:)))
% max(abs(D(:)-D1(:)))
% 
% UH = inituhat(master,mesh1.elcon,UDG,1);
% [Ru1, Rh1, B1, D1, F1, G1, K1, H1] = uequationelemface(master,app1,mesh1,UDG,UH,SDG);
% [Ru, Rh, B, D, F, G, K, H] = uEquationFaceExasim(pde, mesh);
% 
% max(abs(Ru(:)-Ru1(:)))
% max(abs(Rh(:)-Rh1(:)))
% max(abs(B(:)-B1(:)))
% max(abs(D(:)-D1(:)))
% max(abs(F(:)-F1(:)))
% max(abs(G(:)-G1(:)))
% max(abs(K(:)-K1(:)))
% max(abs(H(:)-H1(:)))
% 
% clear app1;
% app1.hybrid = "hdg";
% app1.bcm = [1; 1; 1; 1];
% app1.bcs = [1; 1; 1; 1]*0;
% app1.flg_q = 1;
% app1.flg_p = 1;
% app1.time = 0;
% app1.tdep = 0;
% app1.wave=0;
% app1.fc_u = 1;
% app1.fc_q = 1;
% app1.denseblock=0;
% app1.arg = {};
% app1.flux = "flux";
% app1.source = "source";
% app1.fbou   = "fbou";
% app1.fhat   = "fhat";
% [master,mesh1] = preprocess(master,mesh1,app1.hybrid);
% 
% pde = initializepde("backend");
% pde.porder = 4;             % polynomial degree
% pde.physicsparam = 1;       % unit thermal conductivity
% pde.pgauss = 2*pde.porder;
% pde.nd = 2;
% pde.arg = {};
% pde.bcm = [1; 1; 1; 1];
% pde.bcs = [1; 1; 1; 1]*0;
% pde.denseblock=0;
% pde.flg_p = 0;
% pde.localsolve = 1;
% pde.adjoint = 0;
% pde.flux = "flux";
% pde.source = "source";
% pde.fbou   = "fbou";
% pde.fhat   = "fhat";
% pde.nc = 3;
% 
% master1 = Master(pde);
% mesh1 = hdgmesh(mesh, pde.porder);
% UH = inituhat(master1,mesh1.elcon,UDG,1);
% [UDG,UH]= hdgsolve(master1,mesh1,pde,UDG,UH,[]);
% 
% 
