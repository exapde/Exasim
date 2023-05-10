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
pde.porder = 1;             % polynomial degree
physicsparam = 1:10; % unit thermal conductivity and zero boundary data
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
% run FOM to initialize everything
% [sol,pde,mesh,master,dmd,compilerstr,runstr] = exasim(pde,mesh);
pde = setcompilers(pde);       

% generate input files and store them in datain folder
[pde,mesh,master,dmd] = preprocessing(pde,mesh);

% generate source codes and store them in app folder
gencode(pde);

% compile source codes to build an executable file and store it in app folder
compilerstr = compilecode(pde);

%% Get snapshots
%TODO: would be nice if we could save to multiple out directories 
fileapp = "datain/app.bin";
snapshots = {}
for i = 1:length(physicsparam)
    pde.physicsparam = physicsparam(i);
    writeapp(pde,fileapp,'native');
    runcode(pde);
    sol = fetchsolution(pde,master,dmd);
    snapshots{i} = sol;
end

%% arrange snapshots in a matrix
% only use components related to u, not gradient
n_train = length(physicsparam);
snapshots_mat = zeros(numel(snapshots{1}(:,1,:)), n_train);
for i = 1:n_train
    tmp = snapshots{i}(:,1,:);
    snapshots_mat(:,i) = reshape(tmp, [master.npe*pde.ncu*pde.ne, 1]);
end

%% build POD basis by taking the SVD
% TODO: do correct POD formulation, eigenvecs u of Phi^T Phi and then W = Phi u
% [U,S,V] = svd(snapshots_mat);
[V,D] = eig(snapshots_mat'*snapshots_mat);
[d,ind] = sort(diag(D),'descend');
Ds = D(ind,ind);
Vs = V(:,ind);
Phi_POD = snapshots_mat*Vs;

%% Extract reduced basis. Here we take the 5 leading singular vectors
figure(1); clf; semilogy(diag(Ds),'LineWidth',2); title("Singular values");
p_rb = 3;
% Phi = U(:,1:p_rb);


Phi_RB = snapshots_mat(:,[1, 2, 3]);

%% Initialize Galerkin ROM
% We have Phi...now we should do a quick galerkin test

Phi = Phi_POD(:,1:p_rb);
% Set up runmode and parameter
mu_test = 2.5;
pde.runmode = 0;
pde.physicsparam = mu_test;
writeapp(pde,fileapp,'native')
runcode(pde);
solTrue = fetchsolution(pde,master,dmd);
solTrue = solTrue(:,1,:);


pde.runmode = 20; % run mode > 5: R(u,mu) and J(u,mu) v
writeapp(pde,fileapp,'native')

% Set up newton step initializations
u_rb = zeros(p_rb,1);
mesh.udg = reshape(Phi * u_rb, size(snapshots{1}(:,1,:)));
mesh.dudg = 0*mesh.udg;
writesol(pde,mesh,master,dmd)
runcode(pde);
R0 = getsolution('dataout/out_Ru_test',dmd,master.npe);
Jv0 = getsolution('dataout/out_Jv_test',dmd,master.npe);
b = Phi'*reshape(R0, [master.npe*pde.ncu*pde.ne, 1]); %TODO: probably need to be more careful about sizings here

%% Run Galerkin ROM
% System to solve is Phi^T R(Phi u_rb) = 0
% First step: let's set up 
maxiter = 20;
iter = 1;
JPhi = 0*Phi;
% for k = 1:5
while norm(b) > 1e-8 && iter < maxiter
    Phiu = Phi * u_rb;
    mesh.udg = reshape(Phiu, [master.npe, pde.ncu, pde.ne]);
    for i = 1:p_rb
        % Evaluate JPhi by calling Exasim's Jv against each column of Phi
        mesh.dudg = reshape(snapshots_mat(:,i), [master.npe, pde.ncu, pde.ne]);
        writesol(pde,mesh,master,dmd)
        runcode(pde);
        Jv = reshape(getsolution('dataout/out_Jv_test',dmd,master.npe),[master.npe*pde.ncu*pde.ne, 1]);
        JPhi(:,i) = Jv;
    end
    R = reshape(getsolution('dataout/out_Ru_test',dmd,master.npe),[master.npe*pde.ncu*pde.ne, 1]);
%%%%%%%%% Galerkin
    A =  Phi' * JPhi;
    b = -Phi' * R;
%%%%%%%%%

%%%%%%%%%% Gauss Newton TODO: weighted inner product...
%%% Could do Minv inside code or not...
%     A =  JPhi' * JPhi;
%     b = -JPhi' * R;
%%%%%%%%%%

    du_rb = A\b;
    u_rb = u_rb + du_rb;
    iter = iter+1;
    disp(norm(b));
end
%% Evalute R(u) and J(u)v

% %update app to run just R, Jv
% pde.runmode = 20;
% fileapp = "datain/app.bin";
% writeapp(pde,fileapp,'native');
% 
% % load in u to mesh.udg and v to mesh.dudg
% Uout = getsolution('dataout/out',dmd,master.npe);
% mesh.udg = Uout(:,1:1,:);
% mesh.dudg = mesh.udg;
% % write to binary
% 
% runcode(pde);
% Rout = getsolution('dataout/out_Ru_test',dmd,master.npe);
% Jvout = getsolution('dataout/out_Jv_test',dmd,master.npe);

figure(1); clf; scaplot(mesh, solTrue); title("True solution")
figure(2); clf; scaplot(mesh, Phi*u_rb); title("Galerkin ROM")
% figure(3); scaplot(mesh, Phi*u_rb - solTrue(:)); title("Error")
fprintf('Maximum absolute error: %g\n',max(abs(Phi*u_rb-solTrue(:))));

%%
% visualize the numerical solution of the PDE model using Paraview
% pde.visscalars = {"temperature", 1};                % list of scalar fields for visualization
% pde.visvectors = {"temperature gradient", [2 3 4]}; % list of vector fields for visualization
% dgnodes = vis(sol,pde,mesh);                        % visualize the numerical solution
% x = dgnodes(:,1,:); y = dgnodes(:,2,:); z = dgnodes(:,3,:);
% uexact = sin(pi*x).*sin(pi*y).*sin(pi*z);           % exact solution
% uh = sol(:,1,:);                                    % numerical solution
% fprintf('Maximum absolute error: %g\n',max(abs(uh(:)-uexact(:))));
disp("Done!");
