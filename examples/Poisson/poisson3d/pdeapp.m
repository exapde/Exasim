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
pde.physicsparam = [1 0.0]; % unit thermal conductivity and zero boundary data
pde.tau = 1.0;              % DG stabilization parameter
pde.linearsolvertol = 1e-6; % GMRES tolerance
pde.ppdegree = 1;           % degree of polynomial preconditioner
pde.RBdim = 0;              % reduced basis dimension for preconditioner
pde.GMRESrestart = 20;

% create a grid of 8 by 8 by 8 hexes on the unit cube
[mesh.p,mesh.t] = cubemesh(8,8,8,1);
% mesh.p = mesh2.p';
% mesh.t = mesh2.t';
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(2,:))<1e-8, @(p) abs(p(1,:)-1)<1e-8, @(p) abs(p(2,:)-1)<1e-8, @(p) abs(p(1,:))<1e-8, @(p) abs(p(3,:))<1e-8, @(p) abs(p(3,:)-1)<1e-8};
mesh.boundarycondition = [1;1;1;1;1;1]; % Set boundary condition for each boundary

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh,master,dmd] = exasim(pde,mesh);

% visualize the numerical solution of the PDE model using Paraview
pde.visscalars = {"temperature", 1};                % list of scalar fields for visualization
pde.visvectors = {"temperature gradient", [2 3 4]}; % list of vector fields for visualization
dgnodes = vis(sol,pde,mesh);                        % visualize the numerical solution

%dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);    
x = dgnodes(:,1,:); y = dgnodes(:,2,:); z = dgnodes(:,3,:);
uexact = sin(pi*x).*sin(pi*y).*sin(pi*z);           % exact solution
uh = sol(:,1,:);                                    % numerical solution
fprintf('Maximum absolute error: %g\n',max(abs(uh(:)-uexact(:))));
disp("Done!");

return;

% pde.elemtype = 0;
% master1 = Master(pde);
% mesh1 = hdgmesh(mesh, pde.porder);
% UDG0 = 0*mesh1.dgnodes;
% UDG0(:,4,:) = 0;
% UH0 = inituhat(master1,mesh1.elcon,UDG0,1);
% pde.bcm = [1; 1; 1; 1; 1; 1];
% pde.bcs = [1; 1; 1; 1; 1; 1]*0;
% pde.debugmode=0;
% pde.denseblock=1;
% [UDG,UH]= hdgsolve(master1,mesh1,pde,UDG0,UH0,[]);
% %compareexasim(master1, mesh1, pde);
% 
% x = mesh1.dgnodes(:,1,:); y = mesh1.dgnodes(:,2,:); z = mesh1.dgnodes(:,3,:);
% uexact = sin(pi*x).*sin(pi*y).*sin(pi*z);           % exact solution
% uh = UDG(:,1,:);                                    % numerical solution
% fprintf('Maximum absolute error: %g\n',max(abs(uh(:)-uexact(:))));
% disp("Done!");


