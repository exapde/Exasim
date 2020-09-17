% Specify an Exasim version to run
version = "Version0.1"; 

% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/"  + version + "/Matlab" + "/setpath.m");

%%% Create app struct %%%
app = initializeapp(version);

app.appname = "poisson"; % application name 
app.platform = "cpu";   % choose this option if NVIDIA GPUs are not available 
%app.platform = "gpu";  % choose this option if NVIDIA GPUs are available 
app.mpiprocs = 1;       % number of MPI ranks

%%% Set CPU compiler %%%
%%% use default g++ compiler on your system
app.cpucompiler = "g++";  
%%% or modify the below option to choose a specific CPU compiler   
% app.cpucompiler = "/opt/local/bin/g++-mp-10";

%%% set CPU compiler options %%%
%%% use the below options if blas/lapack libary is NOT in the system search path
% app.cpuflags = "-O2 -ldl -lm -Wl,-rpath, /path/to/blaslapack -L/path/to/blaslapack -lblas -llapack";
%%% use the below options if MKL libary is available on your system 
% app.cpuflags = "-O2 -ldl -lm -Wl,-rpath, /path/to/MKL -L/path/to/MKL -lmkl_intel_lp64 -lmkl_sequential -lmkl_core";

%%% Set MPI compiler if app.mpiprocs > 1 %%%
% app.mpicompiler = "mpicxx"; 
% app.mpirun = "mpirun";
%%% or modify the below options to choose a specific MPI compiler   
% app.mpicompiler = "/opt/local/bin/mpicxx-openmpi-mp"; 
% app.mpirun = "/opt/local/bin/mpirun-openmpi-mp";

%%% Set GPU compiler if app.platform = "gpu" %%%
%%% use default nvcc compiler on your system
% app.gpucompiler = "nvcc";
%%% or modify the below option to choose a specific GPU compiler   
% app.gpucompiler = "/opt/local/bin/nvcc";

%%% set GPU compiler options %%%
%%% use the below options if CUDA libaries are NOT in the system search path
% app.gpuflags = "-Wl,-rpath, /path/to/CUDAlibs -L/path/to/CUDAlibs -lcudart -lcublas";

% Define PDE model, governing equations, and boundary conditions
app.pdemodel = 2;       % (u,q) type
app.porder = 3;         % polynomial degree

% Poisson equation with homogenous Dirichlet condition on a unit square
% q + \nabla u = 0 in \Omega
% \nabla dot flux = source in \Omega
% u = 0 on \partial Omega
% flux = 2*param*q
% source = 2*pi*pi*sin(pi*x)*sin(pi*y);
app.Flux = "flux";      % name of the function defining PDE fluxes
app.Source = "source";  % name of the function defining source term
app.Fbou = "fbou";      % name of the function defining boundary flux
app.Ubou = "ubou";      % name of the function defining boundary value for the solution
app.dt = 0;             % steady-state problem
app.boundaryconditions = [1;1;1;1]; % Set boundary condition for each boundary
app.physicsparam = 1;   % unit thermal conductivity 
app.tau = 1;            % stabilization parameter
app.ncu = 1;            % number of state variables

% mesh struct
m = 8; n = 8;
elemtype = 1;
% create a grid of 8 by 8 on the unit square
[mesh.p,mesh.t] = squaremesh(m,n,1,elemtype);
% expressions for domain boundaries
mesh.bndexpr = {@(p) abs(p(2,:))<1e-8, @(p) abs(p(1,:)-1)<1e-8, @(p) abs(p(2,:)-1)<1e-8, @(p) abs(p(1,:))<1e-8};
% expressions for curved boundaries
mesh.curvedboundary = [0 0 0 0];
mesh.curvedboundaryexpr = {};
% experssions for periodic boundaries
mesh.periodicexpr = {};
% prdexpr = {2, @(p) p(2,:), 4, @(p) p(2,:)};
% prdexpr = {1, @(p) p(1,:), 3, @(p) p(1,:)};
% prdexpr = {1, @(p) p(1,:), 3, @(p) p(1,:); ...
%           2, @(p) p(2,:), 4, @(p) p(2,:)};

% generate curved high-order mesh from the linear mesh
mesh = createhighordermesh(mesh,app);  

% sol struct
sol.UDG = 0*mesh.dgnodes;
sol.UDG(:,3,:) = 0;
sol.ODG = [];
sol.WDG = [];

% run preprocessing to save input data into binary files
[app,master,dmd] = preprocessing(app,mesh,sol);

% generate source codes and store them in app
gencode(app);

% compile source codes to build an executable file and store it in app
compilerstr = compilecode(app);

% run executable file to compute solution and store it in dataout
runstr = runcode(app);

% read solution from dataout and visualize it 
app.vissurf = 1;
sol = vis(sol,app,mesh,master,dmd);

% execute a parametric study by looping over a set of physical parameters
for i = 1:size(paramset,1)
    % set app.physicparam to a new parameter vector
    app.physicparam = paramset(i,:);
    
    % save the modified app struct into the binary file app.bin
    writeapp(app,"datain/app.bin");    
    
    % run executable file to compute solution for the modified app struct
    runcode(app);
    
    % get solution from output files in dataout folder
    Uout = fetchsolution(app,master,dmd);    
    
    % do something with Uout
end



