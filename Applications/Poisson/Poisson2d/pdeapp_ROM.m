% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
ii = ii(end);
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();
pde.cpucompiler="/opt/homebrew/Cellar/llvm@15/15.0.7/bin/clang++";

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";       % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel"; % name of a file defining the PDE model
% pde.enzyme = "ClangEnzyme-15.dylib";

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 2;             % polynomial degree
pde.tau = 1.0;              % DG stabilization parameter
pde.physicsparam=1.0; % placeholder for gencode
% Choose computing platform and set number of processors
%pde.platform = "gpu";      % choose this option if you want to run the C++ code on Nvidia GPUs
pde.mpiprocs = 1;           % number of MPI processors

% create a grid of 8 by 8 on the unit square
[mesh.p,mesh.t] = squaremesh(25,25,1,1);
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
%TODO: would be nice if we could save to multiple out directories  

% Initialize training set
n_train = 10;
parameters = cell(n_train,1);
for itrain = 1:10
    parameters{itrain} = itrain;
end

% Run training cases
snapshots = run_parameter_ensemble(parameters, pde, master, dmd);

%% arrange snapshots in a matrix
Phi = snapshot_struct_to_matrix(snapshots, pde, master);

%% build POD basis by taking the SVD
[Mi, M] = massinv(pde, master,mesh);
weight = M;
[W_POD_full, svals] = snapshot_POD(Phi, weight, pde, master);
figure(1); clf; semilogy(abs(svals)./svals(1),'LineWidth',2); title("Singular values");

%% Extract reduced basis. 
p_rb = 3;
W_RB = Phi(:,[1, 2, 3]);
W_POD = W_POD_full(:,1:p_rb);


W = W_POD; % Take leading POD vectors
% W = W_RB; % RBM basis 
%% Initialize Galerkin ROM
% Set up test set
n_test = 1;
mu_test = cell(n_test,1);
mu_test{1} = 2.5;
snapshots_test = run_parameter_ensemble(mu_test, pde, master, dmd);

%% Run Galerkin ROM
% System to solve is Phi^T R(Phi u_rb) = 0
u_0 = zeros(p_rb,1);
[u_rb, u_rb_history] = run_steady_ROM(mu_test{1}, u_0, W, pde, mesh,master, dmd, "G");


solTrue = snapshots_test{1}(:,1,:);
figure(1); clf; scaplot(mesh, solTrue); title("True solution")
figure(2); clf; scaplot(mesh, W*u_rb); title("ROM")
% figure(3); scaplot(mesh, Phi*u_rb - solTrue(:)); title("Error")
fprintf('Maximum absolute error: %g\n',max(abs(W*u_rb-solTrue(:))));
disp("Done!");
