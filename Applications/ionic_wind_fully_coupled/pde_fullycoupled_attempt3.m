% clear exasim data from memory
clear pde mesh master dmd sol;

% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel_fc_attempt3";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
%pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors

% Physical parameters
Kep = 2e-13;             % mu[1] Recombination coeff - pos and neg ions [m^3/s]
Knp = 2e-13;             % mu[2] Recombination coeff - pos ions and electrons [m^3/s]
mu_p = 2.43e-4;          % mu[3] Pos ion mobility [m^2/(Vs)]
mu_n = 2.7e-4;           % mu[4] Neg mobility [m^2/(Vs)]
De = 0.18;               % mu[5] Electron diffusion coefficient [m^2/s]
Dp = 0.028e-4;           % mu[6] Pos ion diffusion coefficient [m^2/s]
Dn = 0.043e-4;           % mu[7] Neg diffusion coefficient [m^2/s]
Nmax = 1e16;             % mu[8] Max number density for initial charge distribution [particles/m^3]
r0 = 0.0;                % mu[9] r-pos of emitter tip in reference frame [m]
z0 = 0.045;              % mu[10]z-pos of emitter tip in reference frame [m]
s0 = 1e-2;               % mu[11]Std deviation of initial charge distribution [m]
e = 1.6022e-19;          % mu[12]Charge on electron [C]
epsilon = 8.854e-12;     % mu[13]absolute permittivity of air [C^2/(N*m^2)]
Ua = -10e3;              % mu[14]Emitter potential relative to ground [V]
gamma = 0.001;           % mu[15]Secondary electron emission coefficient [1/m]
E_bd = 3e6;              % mu[16]Breakdown E field in air [V/m]
r_tip = 220e-6;          % mu[17] Tip radius of curvature [m]

% Set discretization parameters, physical parameters, and solver parameters
              %     1     2    3     4     5  6   7    8    9   10  11  12   13     14   15     16     17
pde.physicsparam = [Kep, Knp, mu_p, mu_n, De, Dp, Dn, Nmax, r0, z0, s0, e, epsilon, Ua, gamma, E_bd, r_tip];
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
pde.torder = 1;          % time-stepping order of accuracy
pde.nstage = 1;          % time-stepping number of stages
pde.dt = 1.0e-5*ones(1,3);   % time step sizes
pde.visdt = 1.0;        % visualization timestep size
pde.soltime = 1:pde.visdt:length(pde.dt); % steps at which solution are collected
pde.GMRESrestart=200;            % number of GMRES restarts
pde.linearsolveriter=1000;        % number of GMRES iterations
pde.NLiter=3;                   % Newton iterations

[mesh.p,mesh.t] = gmsh2pt(['chen_geom_coarse.msh'],2, 0);

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
ENUM_SYMMETRY_BC = 2;
ENUM_NEEDLE_BC = 3;

mesh.boundaryexpr = {bdry1,
                    bdry2,
                    bdry3,
                    bdry4,
                    bdry5,
                    bdry6,
                    bdry7};

mesh.boundarycondition = [2, 5, 5, 3, 4, 4, 1]; % Set boundary condition for each boundary

load poisson_sol.mat sol;
mesh.vdg = sol;
% call exasim to generate and run C++ code to solve the PDE models
[sol,pde,mesh,master,dmd,compilerstr,runstr] = exasim(pde,mesh);

% % visualize the numerical solution of the PDE model using Paraview
% for m = 1:length(pde)
%     pde{m}.visscalars = {"temperature", 1};  % list of scalar fields for visualization
%     pde{m}.visvectors = {"temperature gradient", [2 3]}; % list of vector fields for visualization
%     pde{m}.visfilename = "dataout" + num2str(m) + "/output";  
%     vis(sol{m},pde{m},mesh{m}); % visualize the numerical solution
% end
% disp("Done!");

for i = 1:size(sol,4)
    sol(:,[2 4 6],:,i) = sol(:,[2 4 6],:,i) + mesh.vdg;
end

% % visualize the numerical solution of the PDE model using Paraview
% pde.visscalars = {"T1", 1, "T2", 2, "T3", 3, "T4", 4};  % list of scalar fields for visualization
pde.visscalars = {"ne", 1, 'phi', 2};  % list of scalar fields for visualization
pde.visvectors = {"grad ne", [3 5], "grad phi", [4 6]}; % list of vector fields for visualization
xdg = vis(sol,pde,mesh); % visualize the numerical solution
disp("Done!");

% mesh.vdg=sol0*0;    % Zero out the solution before passing it into mesh.vdg