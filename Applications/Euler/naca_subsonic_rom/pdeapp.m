% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
ii = ii(end);
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde structure and mesh structure
tmp = load('meshIDG.mat');
meshIDG = tmp.mesh;
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelC";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
%pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.cpucompiler="/opt/homebrew/Cellar/llvm@12/12.0.1_1/bin/clang++";

pde.mpiprocs = 1;              % number of MPI processors

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 2;          % polynomial degree
pde.torder = 1;          % time-stepping order of accuracy
pde.nstage = 1;          % time-stepping number of stages
pde.dt = 0.0003*ones(1,50000);   % time step sizes
pde.saveSolFreq = 1000;          % solution is saved every 10 time steps
pde.soltime = 1000:10000:50000;% steps at which solution are collected
pde.visdt = pde.dt(1);           % visualization timestep size

gam = 1.4;              % specific heat ratio
Minf = 0.5;             % Mach number
rinf = 1.0;             % freestream density
aoa = 3;
alpha = aoa*pi/180;
uinf = cos(alpha);             % freestream horizontal velocity
vinf = sin(alpha);             % freestream vertical velocity
pinf = 1/(gam*Minf^2);  % freestream pressure
rEinf = 0.5+pinf/(gam-1); % freestream energy
pde.physicsparam = [gam Minf rinf uinf vinf rEinf];
pde.tau = 2.0;          % DG stabilization parameter

pde.GMRESrestart=30;  % number of GMRES restarts
pde.linearsolvertol=0.0001; % GMRES tolerance
pde.linearsolveriter=31;  % number of GMRES iterations
pde.precMatrixType=2; % preconditioning type
pde.NLtol = 1e-7;  % Newton tolerance
pde.NLiter = 3;   % number of Newton iterations
% pde.ppdegree = 15;

% read a grid from a file
% [mesh.p,mesh.t] = readmesh('grid.bin',0);
mesh.p = meshIDG.p';
mesh.t = meshIDG.t';
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) sqrt((p(1,:)-.5).^2+p(2,:).^2)<3, @(p) abs(p(1,:))<20};
% mesh.boundarexpr = {@(p) sqrt(p(1,:).^2 + p(2,:).^2)<2,@(p) sqrt(p(1,:).^2 + p(2,:).^2)<20};  
mesh.boundarycondition = [1;2];
% expressions for curved boundaries
mesh.curvedboundary = [1 0];
mesh.curvedboundaryexpr = {@(p) p(2,:).^2-(5*0.01*12*(0.29690*sqrt(abs(p(1,:)))-0.12600*p(1,:)-0.35160*p(1,:).^2+0.28430*p(1,:).^3-0.10150*p(1,:).^4)).^2, @(p) 0};

%%
pde = setcompilers(pde);       
% generate input files and store them in datain folder
[pde,mesh,master,dmd] = preprocessing(pde,mesh);
% generate source codes and store them in app folder
gencode(pde);
% compile source codes to build an executable file and store it in app folder
compilerstr = compilecode(pde);
% call exasim to generate and run C++ code to solve the PDE model
% [sol,pde,mesh,master,dmd,compilerstr,runstr,res] = exasim(pde,mesh);
mesh.porder = pde.porder;
ncu = pde.ncu;
%% Get snapshots
% Define training samples
n_train = 5;
parameters = cell(n_train^2,1);
Ma_range = linspace(0.3, 0.5, n_train);
aoa_range = linspace(0,3, n_train);
mu_ind = 1;
for itrain = 1:n_train
    for jtrain = 1:n_train
        Minf = Ma_range(itrain);
        aoa = aoa_range(jtrain);       
        alpha = aoa*pi/180;
        uinf = cos(alpha);             % freestream horizontal velocity
        vinf = sin(alpha);             % freestream vertical velocity
        pinf = 1/(gam*Minf^2);  % freestream pressure
        rEinf = 0.5+pinf/(gam-1); % freestream energy
        physicsparam = [gam Minf rinf uinf vinf rEinf];
        parameters{mu_ind} = physicsparam;
        mu_ind = mu_ind + 1;
    end
end

% Run training cases
snapshots = run_parameter_ensemble(parameters, pde, master, dmd);

%%
%% arrange snapshots in a matrix
Phi = snapshot_struct_to_matrix(snapshots, pde, master);

%% build POD basis by taking the SVD
% [Mi, M] = massinv(pde, master,mesh);
weight = 1;
[W_POD_full, svals] = snapshot_POD(Phi, weight, pde, master);
figure(1); clf; semilogy(abs(svals)./svals(1),'LineWidth',2); title("Singular values");

%% Extract reduced basis. 
p_rb = 11;
W_RB = Phi(:,[1, 2, 3]);
W_POD = W_POD_full(:,1:p_rb);

%TODO: next step should be to actually access elements for a more accurate
%comparison. 
i_hr = unique(randi(pde.ne, 30, 1));
% ii = sort(ii);
% ii = 1:pde.ne*master.npe;
W = W_POD; % Take leading POD vectors
% W = W_RB; % RBM basis 

% JPhi = r
% JPhi_train = rand([size()]) %what is the best way to evaluate this...;

%%
n_test = 1;
Minf = 0.375;
aoa = 1.15;       
alpha = aoa*pi/180;
uinf = cos(alpha);             % freestream horizontal velocity
vinf = sin(alpha);             % freestream vertical velocity
pinf = 1/(gam*Minf^2);  % freestream pressure
rEinf = 0.5+pinf/(gam-1); % freestream energy
mu_test = cell(n_test,1);
mu_test{1} = [gam Minf rinf uinf vinf rEinf];
snapshots_test = run_parameter_ensemble(mu_test, pde, master, dmd);
solTrue = snapshots_test{1}(:,1:ncu,:,end);

%% Run Galerkin ROM
% System to solve is Ph i^T R(Phi u_rb) = 0

%TODOs: Better initial guess
%
err_G = zeros(p_rb,ncu);
err_G_S = zeros(p_rb,ncu);
err_hr = zeros(p_rb,ncu);
err_L = zeros(p_rb,ncu);
err_L_hr = zeros(p_rb,ncu);
err_0 = zeros(p_rb,ncu);
for p = 1:p_rb

disp("==================================================");
disp("=========RUNNING p" + string(p)  + "==============");
disp("==================================================");
W_tmp = W_POD_full(:,1:p);
% u_0 = zeros(p,1);
% u_0(1) = 1;
u_0 = 1/4 * pinv(W_tmp) * (Phi(:,7) + Phi(:,8) + Phi(:,12) + Phi(:,13));
JW_train = run_parameter_ensemble_RJV(parameters, Phi, W_tmp, pde, master, mesh, dmd);

[u_rb_hr,~]    = run_steady_ROM_hr_reg_test(mu_test{1}, u_0, W_tmp, pde, mesh, master, dmd, "G", i_hr, JW_train);
[u_rb_g, ~]    = run_steady_ROM(mu_test{1}, u_0, W_tmp, pde, mesh,master, dmd, "G",1:(pde.ne));
[u_rb_l, ~]    = run_steady_ROM(mu_test{1}, u_0, W_tmp, pde, mesh,master, dmd, "L",1:(pde.ne));
[u_rb_l_hr, ~] = run_steady_ROM_hr_reg_test(mu_test{1}, u_0, W_tmp, pde, mesh,master, dmd, "L",i_hr,JW_train);
[u_rb_g_S, ~]  = run_steady_ROM(mu_test{1}, u_0, W_tmp, pde, mesh,master, dmd, "G",i_hr);

UDG_ROM_ig    = reshape(W_tmp*u_0, size(solTrue));
UDG_ROM_g    = reshape(W_tmp*u_rb_g, size(solTrue));
UDG_ROM_hr   = reshape(W_tmp*u_rb_hr, size(solTrue));
UDG_ROM_g_S  = reshape(W_tmp*u_rb_g_S, size(solTrue));
UDG_ROM_l    = reshape(W_tmp*u_rb_l, size(solTrue));
UDG_ROM_l_hr = reshape(W_tmp*u_rb_l_hr, size(solTrue));

err_0(p,:)    = computeerror_ROM(solTrue, UDG_ROM_ig,mesh,master);
err_G(p,:)    = computeerror_ROM(solTrue, UDG_ROM_g,mesh,master);
err_G_S(p,:)  = computeerror_ROM(solTrue, UDG_ROM_g_S,mesh,master);
err_L(p,:)    = computeerror_ROM(solTrue, UDG_ROM_l,mesh,master);
err_hr(p,:)   = computeerror_ROM(solTrue, UDG_ROM_hr,mesh,master);
err_L_hr(p,:) = computeerror_ROM(solTrue, UDG_ROM_l_hr,mesh,master);
end


%% Postprocessing
figure(5); clf; hold on; 
semilogy(err_G(1:end,end),'LineWidth',2);
semilogy(err_L(1:end,end),'LineWidth',2);
semilogy(err_G_S(1:end,end),'--','LineWidth',2);
semilogy(err_hr(1:end,end),'--o','LineWidth',2);
semilogy(err_L_hr(1:end,end),'--o','LineWidth',2);
legend(["Galerkin", "LSPG", "Reduced Galerkin", "Galerkin Expand Jacobian", "LSPG Expand Jacobian"])

%%
% visualize the numerical solution of the PDE model using Paraview
pde.visscalars = {"density", 1, "energy", 4};  % list of scalar fields for visualization
pde.visvectors = {"momentum", [2, 3]}; % list of vector fields for visualization
xdg = vis(sol,pde,mesh); % visualize the numerical solution
disp("Done!");
