% clear exasim data from memory
clear pde mesh master dmd sol;

% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel_full_model";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
%pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors

% Physical parameters
r0 = 0.0;                % mu[9] r-pos of emitter tip in reference frame [m]
z0 = 0.0;                % mu[10]z-pos of emitter tip in reference frame [m]
s0 = 25e-6;               % mu[11]Std deviation of initial charge distribution [m]
Nmax = 1e16;             % mu[8] Max number density for initial charge distribution [particles/m^3]
e = 1.6022e-19;          % mu[12]Charge on electron [C]
epsilon0 = 8.854e-12;     % mu[13]absolute permittivity of air [C^2/(N*m^2)] or [m-3 kg-1 s4 A2]
Ua = -10e3;              % mu[14]Emitter potential relative to ground [V]
gamma = 0.001;           % mu[15]Secondary electron emission coefficient [1/m]
E_bd = 3e6;              % mu[16]Breakdown E field in air [V/m]
r_tip = 220e-6;          % mu[17] Tip radius of curvature [m]
n_ref = epsilon0*E_bd/(e*r_tip);     % Yes, this is redundant and could be recomputed from the above variables. But it saves having to recompute it each time in the functions.

P = 101325; % Pa
V = 1; % m^3
T = 273.15; % K
k_b = 1.380649e-23; % m2 kg s-2 K-1
N = P*V/(k_b*T);     % Neutral number density

% Set discretization parameters, physical parameters, and solver parameters
              %     1   2    3   4    5      6     7     8     9      10    11    12
pde.physicsparam = [r0, z0, s0, Nmax, e, epsilon0, Ua, gamma, E_bd, r_tip, n_ref, N];
pde.tau = 1.0;           % DG stabilization parameter

% set indices to obtain v from the solutions of the other PDE models 
% first column : model index
% second column: solution index
pde.porder = 3;          % polynomial degree
% pde.vindx = [1 1; 1 2; 1 3]; % first column -> model, second column -> solution index for that model. Originally this used to be "v"
% 4, 5 6, 7 would be grad ne_x, grad ne_y, grad ne_x, etc
% pde.subproblem = 1;

pde.NLtol = 1.0e-14;
pde.linearsolvertol = 1.0e-14;
pde.ppdegree = 0;      % polynomial preconditioner degree -> set to 0 because we aren't using the PPC
% pde.precMatrixType = 2;

% solver parameters
% Notes on computing the nondimensional timestep:
% - mue at Ebd and N computed at STP = 2e-4

mue_ref = 2e-4;
dt_sec = 1.0e-5;     % This is the dimensional timestep, in seconds

% dt_star = dt_sec*mue_ref*E_bd/r_tip;     % This is the nondimensional timestep, delta t*
dt_star = 1e-3;     % Setting it to 1e-3 for testing
% disp('The nondimensional timestep is')
% disp(dt_star)
pde.torder = 1;          % time-stepping order of accuracy
pde.nstage = 1;          % time-stepping number of stages
pde.dt = dt_star*ones(1,3);   % time step sizes
pde.visdt = 1.0;        % visualization timestep size
pde.soltime = 1:pde.visdt:length(pde.dt); % steps at which solution are collected
pde.GMRESrestart=50;            % number of GMRES restarts
pde.linearsolveriter=1000;        % number of GMRES iterations
pde.NLiter=3;                   % Newton iterations

[mesh.p,mesh.t] = gmsh2pt(['chen_geom_coarse.msh'],2, 0);

% expressions for domain boundaries
eps = 1e-4;
xmin = min(mesh.p(1,:));
xmax = max(mesh.p(1,:));
ymin = min(mesh.p(2,:));
ymax = max(mesh.p(2,:));
x2 = 0.017;
x3 = 0.015;

bdry1 = @(p) (p(1,:) < xmin+eps);    % Chen bdry 2 - symmetry boundary
bdry2 = @(p) (p(1,:) > xmax - eps);  % Chen bdry 5 - +r farfield
bdry3 = @(p) (p(2,:) > ymax - eps);  % Chen bdry 6 - +z farfield
bdry4 = @(p) (p(2,:) < ymin+eps) && (p(1,:) < x3+eps);   % Chen bdry 3 - outflow boundary
bdry5 = @(p) (p(2,:) < ymin+eps) && (p(1,:) > x2-eps);   % grounded boundary - Chen bdry 4 - ground plane
bdry6 = @(p) (p(1,:) < .02) && (p(1,:) > .01);           % grounded boundary - Chen bdry 4 - cylinder
bdry7 = @(p) (p(1,:) < x2+eps);                          % Chen bdry 1 - needle

mesh.boundaryexpr = {bdry1,
                    bdry2,
                    bdry3,
                    bdry4,
                    bdry5,
                    bdry6,
                    bdry7};

mesh.boundarycondition = [2, 5, 5, 3, 4, 4, 1]; % Set boundary condition for each boundary

% load poisson_sol.mat sol;
% mesh.vdg = sol;
% call exasim to generate and run C++ code to solve the PDE models
[sol,pde,mesh,master,dmd,compilerstr,runstr] = exasim(pde,mesh);

% for i = 1:size(sol,4) % Last field is the number of timesteps
%     sol(:,[4 8 12],:,i) = sol(:,[4 8 12],:,i) + mesh.vdg;
% end

% visualize the numerical solution of the PDE model using Paraview
pde.visscalars = {"ne", 1, "np", 2, "nn", 3, 'phi', 4};  % list of scalar fields for visualization
pde.visvectors = {"grad ne", [5 9], "grad np", [6 10], "grad nn", [7 11], "grad phi", [8 12]}; % list of vector fields for visualization
xdg = vis(sol,pde,mesh); % visualize the numerical solution
disp("Done!");
