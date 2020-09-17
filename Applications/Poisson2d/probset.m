% specify an Exasim version to run
version = 'Version0.1'; 

% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
tmp = strfind(cdir(ii(end):end), "/");
up = length(tmp);
codedir = "";
for i = 1:up
    codedir = codedir + "../";
end
versiondir = codedir  + version;
run(versiondir + "/setup.m");

% mesh struct
ngrid = 9;
elemtype = 1;
% create a grid
[mesh.p,mesh.t] = squaremesh(ngrid,ngrid,1,elemtype);
% expressions for domain boundaries
mesh.bndexpr = {@(p) abs(p(2,:))<1e-8, @(p) abs(p(1,:)-1)<1e-8, @(p) abs(p(2,:)-1)<1e-8, @(p) abs(p(1,:))<1e-8};
% experssions for periodic boundaries
mesh.periodicexpr = [];
% prdexpr = {2, @(p) p(2,:), 4, @(p) p(2,:)};
% prdexpr = {1, @(p) p(1,:), 3, @(p) p(1,:)};
% prdexpr = {1, @(p) p(1,:), 3, @(p) p(1,:); ...
%           2, @(p) p(2,:), 4, @(p) p(2,:)};

% app struct
clear app;
app = initializeapp(version);

app.appname = "poi";   % application name 
app.platform = "cpu";  % choose this option if NVIDIA GPUs are not available 
%app.platform = "gpu";  % choose this option if NVIDIA GPUs are available 
app.mpiprocs = 1;       % number of MPI ranks

% Define PDE model, governing equations, and boundary conditions
app.pdemodel = 2;       % (u,q) type
app.porder = 3;         % polynomial degree
app.Flux = "flux";      % name of the script defining PDE fluxes
app.Source = "source";  % name of the script defining source term
app.Fbou = "fbou";      % name of the script defining boundary flux
app.Ubou = "ubou";      % name of the script defining boundary value for the solution
app.dt = 0;             % steady-state problem
app.boundaryconditions = [1;1;1;1]; % Set boundary condition for each boundary
app.physicsparam = 1;   % unit thermal conductivity 
app.tau = 1;            % stabilization parameter

% sol struct
sol.dgnodes = createnodes(mesh.p,mesh.t,app.porder);
sol.UDG = 0*sol.dgnodes;
sol.UDG(:,3,:) = 0;
sol.ODG = [];
sol.WDG = [];

% run preprocessing to save input data into binary files
[app,master,dmd] = preprocessing(app,mesh,sol);

% generate source codes 
gencode(app);

% compile source codes to build an executable application
compilerstr = compilecode(app);

% run executable to compute solution and store it in dataout
runcode(app);

% read solution from dataout and visualize it 
sol = vis(sol,app,mesh,master,dmd);


