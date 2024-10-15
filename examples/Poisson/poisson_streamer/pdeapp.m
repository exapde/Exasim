% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";       % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel"; % name of a file defining the PDE model

% Choose computing platform and set number of processors
pde.mpiprocs = 1;           % number of MPI processors=
pde.hybrid = 1;             % 0 -> LDG, 1-> HDG
pde.debugmode = 0;

<<<<<<< HEAD
pde.platform = "gpu";
=======
pde.platform = "cpu";
>>>>>>> origin/master
pde.GMRESortho = 0;
pde.ppdegree = 20;

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 2;             % polynomial degree
pde.pgauss = 2*pde.porder;  % gauss quad order
% pde.physicsparam = [1 0.0]; % unit thermal conductivity and zero boundary data
pde.tau = 1.0;              % DG stabilization parameter
pde.linearsolvertol = 1e-6; % GMRES tolerance
pde.NLiter = 1;

pde.RBdim = 0;              % reduced basis dimension for preconditioner
pde.linearsolveriter = 1000;
pde.GMRESrestart = 500;

kappa = 1;
l_ref = 1e-4;   % m
mu_ref = 0.05;  % m^2/(V-s)
E_ref = 3e6;    % V/m
e_eps0 = 1.80955e-8;    % Quantity e/epsilon0
phi0 = 18.75e3;     % V
N0 = 5e18; % 1/m^3
z0 = 1e-2; % m
sigma0 = 4e-4;  % m
n_background = 1e13;    % 1/m^3

%          1       2      3       4      5   6    7     8       9
pde.physicsparam = [kappa, l_ref, mu_ref, E_ref, e_eps0, phi0, N0, z0, sigma0];

% mesh = mkmesh_streamer_gmsh(porder, "streamer_380k.msh");
% mesh0 = mkmesh_streamer_gmsh(pde.porder, "streamer_16k_fixed.msh");

<<<<<<< HEAD
[p,t] = gmshcall("streamer_61k.msh", 2, 0);
=======
[p,t] = gmshcall("streamer_16k_fixed.msh", 2, 0);
>>>>>>> origin/master

% Normalization
xmax = max(p(1,:));
p = p/xmax * 125;

mesh.p = p;
mesh.t = t;
% mesh.f = mesh.f';
% mesh.t2f = mesh.t2f';
% mesh.fcurved = mesh.fcurved';
% mesh.tcurved = mesh.tcurved';

% mesh.boundaryexpr = [1;2;3;4]; % Set boundary condition for each boundary
% mesh.boundaryexpr = mesh0.bndexpr;
mesh.boundaryexpr = {@(p) abs(p(2,:))<1e-8, @(p) abs(p(1,:)-125)<1e-8, @(p) abs(p(2,:)-125)<1e-8, @(p) abs(p(1,:))<1e-8};
%mesh.periodicexpr = [0 0 0 0];
mesh.boundarycondition = [1;2;3;4]; % Set boundary condition for each boundary

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh,master,dmd] = exasim(pde,mesh);
mesh.porder = pde.porder;
mesh.dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);
<<<<<<< HEAD
% figure(1); clf; scaplot(mesh,sol(:,1,:),[],2); axis on; axis equal; axis tight;
=======
figure(1); clf; scaplot(mesh,sol(:,1,:),[],2); axis on; axis equal; axis tight;
>>>>>>> origin/master

% pde.platform = "cpu";
% pde.GMRESortho = 0;
% pde.ppdegree = 0;
% [sol,pde,mesh] = exasim(pde,mesh);
% return;
<<<<<<< HEAD

save 'poissonIC61k.mat' sol

normE = sqrt(sol(:,2,:).^2 + sol(:,3,:).^2);
figure(); scaplot(mesh,normE,[],0,0); axis equal; axis tight; colormap jet;
=======
>>>>>>> origin/master
