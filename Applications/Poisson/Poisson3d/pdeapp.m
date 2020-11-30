% Specify an Exasim version to run
version = "Version0.1";

% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim(version);

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";       % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel"; % name of a file defining the PDE model

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;             % polynomial degree
pde.physicsparam = [1 0.0]; % unit thermal conductivity and zero boundary data
pde.tau = 1.0;              % DG stabilization parameter

% Choose computing platform and set number of processors
%pde.platform = "gpu";      % choose this option if you want to run the C++ code on Nvidia GPUs
pde.mpiprocs = 2;           % number of MPI processors

% create a grid of 8 by 8 by 8 hexes on the unit cube
[mesh.p,mesh.t] = cubemesh(8,8,8,0);
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(2,:))<1e-8, @(p) abs(p(1,:)-1)<1e-8, @(p) abs(p(2,:)-1)<1e-8, @(p) abs(p(1,:))<1e-8, @(p) abs(p(3,:))<1e-8, @(p) abs(p(3,:)-1)<1e-8};
mesh.boundarycondition = [1;1;1;1;1;1]; % Set boundary condition for each boundary

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh,master,dmd,compilerstr,runstr] = exasim(pde,mesh);

% visualize the numerical solution of the PDE model using Paraview
pde.visscalars = {"temperature", 1};                % list of scalar fields for visualization
pde.visvectors = {"temperature gradient", [2 3 4]}; % list of vector fields for visualization
dgnodes = vis(sol,pde,mesh);                        % visualize the numerical solution
x = dgnodes(:,1,:); y = dgnodes(:,2,:); z = dgnodes(:,3,:);
uexact = sin(pi*x).*sin(pi*y).*sin(pi*z);           % exact solution
uh = sol(:,1,:);                                    % numerical solution
fprintf('Maximum absolute error: %g\n',max(abs(uh(:)-uexact(:))));
disp("Done!");
