% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel_poisson";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
pde.platform = "cpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;             % number of MPI processors
pde.hybrid = 1;               % 0 -> LDG, 1 -> HDG
pde.debugmode = 0;
pde.nd = 2;
pde.elemtype = 1;
pde.nodetype = 1;

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;             % polynomial degree
pde.pgauss = 2*pde.porder;
pde.tau = 1.0;              % DG stabilization parameter
pde.linearsolvertol = 1e-8; % GMRES tolerance
pde.preconditioner = 1;
pde.GMRESrestart = 50;
pde.ppdegree = 1;          % degree of polynomial preconditioner
pde.RBdim = 0;

% create a grid of 8 by 8 on the unit square
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
pde.physicsparam = [1 theta a1 a2 a];  

% L_i = int phi_i dx
% rho ->  theta = int rho dx / |Omega| 
% F = rho/theta
% int F dx =  int rho/theta dx = int rho dx/ int rho dx / |Omega| = |Omega|

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh,master,dmd,comstr,runstr] = exasim(pde,mesh);

% Apply to the cone problem
% 1. Calculate rho from the flow field 
% 2. Calculate theta
% 3. calculate F = rho/theta 
% 4. calculate grad rho = (drhodx, drhody) 
% 5. mesh.vdg(:,1,:) = F, mesh.vdg(:,2,:) = drhodx/theta, mesh.vdg(:,3,:) = drhody/theta
% 6. Solve the Poisson equation to obtain sol = (u, dudx, dudy)
% 7. mesh.vdg(:,4,:) = sol(:,2,:), mesh.vdg(:,5,:) = sol(:,3,:)
% 8. Solve the transport equation with IC u = x 
% 8. mesh.udg = mesh.dgnodes(:,1,:) and mesh1.dgnodes(:,1,:) = solt(:,1,:,end)
% 9. Solve the transport equation with IC u = y 
% 9. mesh.udg = mesh.dgnodes(:,2,:) and mesh1.dgnodes(:,2,:) = solt(:,1,:,end)
