% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Install/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
%pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors
pde.hybrid = 1;                % 0 -> LDG, 1 -> HDG
pde.debugmode = 1;

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;          % polynomial degree
pde.physicsparam = [1 1];    % unit thermal conductivity
pde.tau = 1.0;           % DG stabilization parameter
pde.linearsolvertol = 1e-3; % GMRES tolerance
pde.GMRESrestart = 100;

% create a grid of 8 by 8 on the unit square
[mesh.p,mesh.t] = squaremesh(4,4,1,1);
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(2,:))<1e-8, @(p) abs(p(1,:)-1)<1e-8, @(p) abs(p(2,:)-1)<1e-8, @(p) abs(p(1,:))<1e-8};
mesh.boundarycondition = [2;1;2;2]; % Set boundary condition for each boundary

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh] = exasim(pde,mesh);

% % visualize the numerical solution of the PDE model using Paraview
% pde.visscalars = {"temperature", 1};  % list of scalar fields for visualization
% pde.visvectors = {"temperature gradient", [2 3]}; % list of vector fields for visualization
% xdg = vis(sol,pde,mesh); % visualize the numerical solution
% disp("Done!");

mesh.porder = pde.porder;
mesh.dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);
figure(1); clf; scaplot(mesh,sol(:,1,:),[],2,1); axis on; axis equal; axis tight;
        
pde.elemtype = 1;
master1 = Master(pde);
mesh1 = hdgmesh(mesh, pde.porder);
UDG0 = 0*mesh1.dgnodes;
UDG0(:,1,:) = 1;
UDG0(:,3,:) = 0;
UH0 = inituhat(master1,mesh1.elcon,UDG0,1);
pde.bcm = [2; 1; 2; 2];
pde.bcs = [1; 1; 1; 1];
[UDG,UH]= hdgsolve(master1,mesh1,pde,UDG0,UH0,[]);
%compareexasim(master1, mesh1, pde);


