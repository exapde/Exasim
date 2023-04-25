% clear exasim data from memory
clear pde mesh master dmd sol;

% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel_poisson";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
%pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors
pde.tau = 1.0;           % DG stabilization parameter

pde.porder = 3;          % polynomial degree
pde.NLtol = 1.0e-6;
pde.linearsolvertol = 1.0e-8;
pde.ppdegree = 0;      % polynomial preconditioner degree -> set to 0 because we aren't using the PPC

[mesh.p,mesh.t] = gmshcall(pde, "chen_geom_coarse2.msh", 2, 0);

% expressions for domain boundaries
eps = 1e-4;
xmin = min(mesh.p(1,:));
xmax = max(mesh.p(1,:));
ymin = min(mesh.p(2,:));
ymax = max(mesh.p(2,:));
% x1 = 4.5e-3;
x2 = 0.017;
x3 = 0.015;

bdry1 = @(p) (p(1,:) < xmin+eps);    % axis symmetric boundary            
bdry2 = @(p) (p(1,:) > xmax - eps);  % open boundary 1                                    
bdry3 = @(p) (p(2,:) > ymax - eps);  % open boundary 2                                     
bdry4 = @(p) (p(2,:) < ymin+eps) && (p(1,:) < x3+eps);   % grounded boundary - open
bdry5 = @(p) (p(2,:) < ymin+eps) && (p(1,:) > x2-eps);   % grounded boundary                                    
bdry6 = @(p) (p(1,:) < .02) && (p(1,:) > .01);                            % grounded boundary - cylinder
bdry7 = @(p) (p(1,:) < x2+eps);                          % needle tip          

ENUM_GROUND_BC = 1;
ENUM_NO_FLUX_BC = 2;
ENUM_NEEDLE_BC = 3;

mesh.boundaryexpr = {bdry1,
                    bdry2,
                    bdry3,
                    bdry4,
                    bdry5,
                    bdry6,
                    bdry7};

mesh.boundarycondition = [ENUM_NO_FLUX_BC,...
                        ENUM_NO_FLUX_BC,...
                        ENUM_NO_FLUX_BC,...
                        ENUM_GROUND_BC,...
                        ENUM_GROUND_BC,...
                        ENUM_GROUND_BC,...
                        ENUM_NEEDLE_BC]; % Set boundary condition for each boundary

% call exasim to generate and run C++ code to solve the PDE models
[sol,pde,mesh,master,dmd,compilerstr,runstr] = exasim(pde,mesh);

isol = readsolstruct('datain/sol.bin');
dgnodes = reshape(isol.xdg, [size(mesh.xpe,1) 2 size(mesh.t,2)]);
save dgnodes.mat dgnodes;
save sol_poi.mat sol;

% visualize the numerical solution of the PDE model using Paraview
% pde.visscalars = {"T1", 1, "T2", 2, "T3", 3, "T4", 4};  % list of scalar fields for visualization
pde.visscalars = {"T1", 1};  % list of scalar fields for visualization
pde.visvectors = {"temperature gradient", [2 3]}; % list of vector fields for visualization
xdg = vis(sol,pde,mesh); % visualize the numerical solution
disp("Done!");