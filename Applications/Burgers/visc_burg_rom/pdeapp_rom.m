% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
ii = ii(end);
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
%pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;          % polynomial degree
kappa = 0.1;
pde.physicsparam = [kappa];    % unit thermal conductivity
pde.tau = 1.0;           % DG stabilization parameter
pde.precMatrixType=2;
% create a grid of 8 by 8 on the unit square
[mesh.p,mesh.t] = squaremesh(64,64,1,1);
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(2,:))<1e-8, @(p) abs(p(1,:)-1)<1e-8, @(p) abs(p(2,:)-1)<1e-8, @(p) abs(p(1,:))<1e-8};
mesh.boundarycondition = [1;1;1;1]; % Set boundary condition for each boundary
mesh.porder = pde.porder;

% call exasim to generate and run C++ code to solve the PDE model
% run FOM to initialize structures
pde = setcompilers(pde);       

% generate input files and store them in datain folder
[pde,mesh,master,dmd] = preprocessing(pde,mesh);

% generate source codes and store them in app folder
gencode(pde);

% compile source codes to build an executable file and store it in app folder
compilerstr = compilecode(pde);

%% Get snapshots

% Define training samples
n_train = 5;
parameters = cell(n_train,1);
tmp_mu = linspace(0, 0.25, n_train+1);
for itrain = 1:n_train
    parameters{itrain} = tmp_mu(itrain+1);
end

% Run training cases
snapshots = run_parameter_ensemble(parameters, pde, master, dmd);
%% arrange snapshots in a matrix
Phi = snapshot_struct_to_matrix(snapshots, pde, master);

%% build POD basis by taking the SVD
[Mi, M] = massinv(pde, master,mesh);
weight = 1;
[W_POD_full, svals] = snapshot_POD(Phi, weight, pde, master);
figure(1); clf; semilogy(abs(svals)./svals(1),'LineWidth',2); title("Singular values");


%% Extract reduced basis. 
p_rb = 5;
W_RB = Phi(:,[1, 2, 3]);
W_POD = W_POD_full(:,1:p_rb);

%TODO: next step should be to actually access elements for a more accurate
%comparison. 
i_hr = unique(randi(pde.ne, 10, 1));
% ii = sort(ii);
% ii = 1:pde.ne*master.npe;
W = W_POD; % Take leading POD vectors
% W = W_RB; % RBM basis 

% JPhi = r
% JPhi_train = rand([size()]) %what is the best way to evaluate this...;

%% Initialize Galerkin ROM
% Set up test set
n_test = 1;
mu_test = cell(n_test,1);
mu_test{1} = 0.213;
snapshots_test = run_parameter_ensemble(mu_test, pde, master, dmd);
solTrue = snapshots_test{1}(:,1,:);

%% Run Galerkin ROM
% System to solve is Ph i^T R(Phi u_rb) = 0
err_G = zeros(p_rb,1);
err_G_S = zeros(p_rb,1);
err_hr = zeros(p_rb,1);
err_L = zeros(p_rb,1);
for p = 1:p_rb
disp("==================================================");
disp("=========RUNNING p" + string(p)  + "==============");
disp("==================================================");
W_tmp = W_POD_full(:,1:p);
u_0 = zeros(p,1);
JW_train = run_parameter_ensemble_RJV(parameters, Phi, W_tmp, pde, master, mesh, dmd);

[u_rb_hr,~] = run_steady_ROM_hr_reg_test(mu_test{1}, u_0, W_tmp, pde, mesh, master, dmd, "G", i_hr, JW_train);
[u_rb_g, ~] = run_steady_ROM(mu_test{1}, u_0, W_tmp, pde, mesh,master, dmd, "G",1:(pde.ne));
% [u_rb_l, ~] = run_steady_ROM(mu_test{1}, u_0, W_tmp, pde, mesh,master, dmd, "L");
[u_rb_g_S, ~] = run_steady_ROM(mu_test{1}, u_0, W_tmp, pde, mesh,master, dmd, "G",i_hr);

UDG_ROM_g = reshape(W_tmp*u_rb_g, size(solTrue));
UDG_ROM_hr = reshape(W_tmp*u_rb_hr, size(solTrue));
UDG_ROM_g_S = reshape(W_tmp*u_rb_g_S, size(solTrue));
% UDG_ROM_l = reshape(W_tmp*u_rb_l, size(solTrue));
err_G(p) = computeerror_ROM(solTrue, UDG_ROM_g,mesh,master);
err_G_S(p) = computeerror_ROM(solTrue, UDG_ROM_g_S,mesh,master);
% err_L(p) = computeerror_ROM(solTrue, UDG_ROM_l,mesh,master);
err_hr(p) = computeerror_ROM(solTrue, UDG_ROM_hr,mesh,master);
end

%% Postprocessing
solTrue = snapshots_test{1}(:,1,:);
figure(1); clf; scaplot(mesh, solTrue); title("True solution")
% figure(2); clf; scaplot(mesh, W_tmp*u_rb_l); title("ROM")
figure(3); clf; scaplot(mesh, abs(W_tmp*u_rb_g - solTrue(:))); title("|u_{Galerkin} - u_h|")
% figure(4); clf; scaplot(mesh, abs(W_tmp*u_rb_l - solTrue(:))); title("|u_{lspg} - u_h|")
%%
figure(5); clf; hold on; 
semilogy(err_G,'LineWidth',2);
% semilogy(err_L,'LineWidth',2);
semilogy(err_G_S,'-','LineWidth',2);
semilogy(err_hr,'--o','LineWidth',2);

grid on;
title("ROM error"); 
% legend(["Galerkin", "LSPG"]);
set(gca,'YScale','log')
xlabel("ROM Dimension")
% UDG_ROM = reshape(W*u_rb, size(solTrue));

% computeerror_ROM(solTrue, UDG_ROM,mesh,master)

% % call exasim to generate and run C++ code to solve the PDE model
% [sol,pde,mesh,master,dmd,compilerstr,runstr,res] = exasim(pde,mesh);
% mesh.porder = pde.porder
% 
% x = mesh.dgnodes(:,1,:); y = mesh.dgnodes(:,2,:);
% uexact = x.*y.*tanh((1-x)/kappa).*tanh((1-y)/kappa);  
% 
% % visualize the numerical solution of the PDE model using Paraview
% pde.visscalars = {"temperature", 1};  % list of scalar fields for visualization
% pde.visvectors = {"temperature gradient", [2 3]}; % list of vector fields for visualization
% xdg = vis(sol,pde,mesh); % visualize the numerical solution
% disp("Done!");
% 
%         