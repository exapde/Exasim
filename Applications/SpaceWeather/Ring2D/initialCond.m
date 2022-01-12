% clear all
% close all
% clc

% Specify an Exasim version to run
version = "Version0.1";
setenv('LD_LIBRARY_PATH', ':/usr/bin');

% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim(version);

% Choose computing platform and set number of processors
% pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors

[pde,mesh] = pdeparams(pde,mesh);

% generate input files and store them in datain folder
[pde,mesh,master,dmd] = preprocessing(pde,mesh);

% sol = mesh.udg;


nElem = size(mesh.dgnodes,3);
nPoints = size(mesh.dgnodes,1);
nDim = 4;
solIC = zeros([nPoints,nDim,nElem]);
% sol2 = zeros([nPoints,nDim,nElem]);

% solIC = initConditionSteady(mesh,solIC,pde.physicsparam);
solIC = initialCondition(mesh,solIC,pde.physicsparam);

% for iElem = 1:nElem
%         for iPoint = 1:nPoints
%             xg = mesh.dgnodes(iPoint,:,iElem);
%             sol(iPoint,:,iElem) = initUfunction(xg,pde.physicsparam,pde.externalparam);
%         end
% end

mesh.nd = master.dim;
mesh.plocal = master.xpe;
mesh.tlocal = master.telem;
mesh.porder = pde.porder;


% visualize the numerical solution of the PDE model using Paraview
pde.visscalars = {"density", 1, "Temperature", 4};  % list of scalar fields for visualization
pde.visvectors = {"velocity", [2, 3]}; % list of vector fields for visualization
xdg = vis(solIC,pde,mesh); % visualize the numerical solution
disp("Done!");
