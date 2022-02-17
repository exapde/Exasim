% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Choose computing platform and set number of processors
% pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model
[pde,mesh] = pdeparams(pde,mesh);

% search compilers and set options
pde = setcompilers(pde);       

% generate input files and store them in datain folder
[pde,mesh,master,dmd] = preprocessing(pde,mesh);
% save('inputFiles','pde','mesh','master','dmd','-v7.3');

% generate source codes and store them in app folder
gencode(pde);

% compile source codes to build an executable file and store it in app folder
compilerstr = compilecode(pde);

% run executable file to compute solution and store it in dataout folder
runstr = runcode(pde);

% get solution from output files in dataout folder
sol = fetchsolution(pde,master,dmd);
sol(:,2,:,:) = sol(:,2,:,:)./sqrt(exp(sol(:,1,:,:)));
sol(:,3,:,:) = sol(:,3,:,:)./sqrt(exp(sol(:,1,:,:)));
sol(:,4,:,:) = sol(:,4,:,:)./sqrt(exp(sol(:,1,:,:)));

% % visualize the numerical solution of the PDE model using Paraview
pde.visscalars = {"density", 1, "temperature", 4};  % list of scalar fields for visualization
pde.visvectors = {"velocity", [2, 3]}; % list of vector fields for visualization
xdg = vis(sol,pde,mesh); % visualize the numerical solution
disp("Done!");

disp("Done!");