% Specify an Exasim version to run
version = "Version0.1";

% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
versiondir = cdir(1:(ii+5)) + "/"  + version + "/Matlab";
run(versiondir + "/setpath.m");

% Define a PDE model: governing equations, initial solutions, and boundary conditions
app = initializeapp(version);  % initialize app struct
app.appname = "poisson";       % application name
app.pdemodel = "ModelD";       % ModelC, ModelD, ModelW
app.pdemodelfile = "pdemodel"; % name of a file defining the PDE model

% Set discretization parameters, physical parameters, and solver parameters
app.porder = 3;          % polynomial degree
app.dt = 0;              % steady-state problem
app.boundaryconditions = [1;1;1;1]; % Set boundary condition for each boundary
app.physicsparam = 1;    % unit thermal conductivity
app.tau = 1.0;           % DG stabilization parameter

% Set compilers and options
app.platform = "cpu";    % choose this option if NVIDIA GPUs are not available
%app.platform = "gpu";   % choose this option if NVIDIA GPUs are available
app.mpiprocs = 1;        % number of MPI processors
app.cpucompiler = "g++"; % Clang/GNU C++ compiler
if app.mpiprocs>1        % set MPI compiler if app.mpiprocs > 1
    app.mpicompiler = "/opt/local/bin/mpicxx-openmpi-mp"; % path to a MPI compiler
    app.mpirun = "/opt/local/bin/mpirun-openmpi-mp";
end
if app.platform == "gpu"
    app.gpucompiler = "nvcc"; % Nvidia CUDA compiler (AMD GPU will be supported in the next version)
end
checkcompilers(app);     % check if the compilers are available

% generate a linear mesh
mesh = initializemesh();
m = 8; n = 8; elemtype = 1;
% create a grid of 8 by 8 on the unit square
[mesh.p,mesh.t] = squaremesh(m,n,1,elemtype);
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(2,:))<1e-8, @(p) abs(p(1,:)-1)<1e-8, @(p) abs(p(2,:)-1)<1e-8, @(p) abs(p(1,:))<1e-8};
% expressions for curved boundaries
mesh.curvedboundary = [0 0 0 0];
mesh.curvedboundaryexpr = {};
% experssions for periodic boundaries
mesh.periodicexpr = {};

% generate input files and store them in datain folder
[app,mesh,master,dmd] = preprocessing(app,mesh);

% generate source codes and store them in app folder
gencode(app);

% compile source codes to build an executable file and store it in app folder
compilerstr = compilecode(app);

% run executable file to compute solution and store it in dataout folder
runstr = runcode(app);

% get solution from output files in dataout folder
sol = fetchsolution(app,master,dmd);

% perform visualization
app.paraview = "/Applications/ParaView-5.8.1.app/Contents/MacOS/paraview";
app.visfilename = "dataout/output";
app.visscalars = {"temperature", 1};
app.visvectors = {"temperature gradient", [2 3]};
app.viselem = 1:size(mesh.t,2);
vis(sol,app,mesh,master); % visualize the numerical solution
disp("Done!");

