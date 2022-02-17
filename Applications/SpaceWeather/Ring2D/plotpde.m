close all

% Specify an Exasim version to run
version = "Version0.1";

% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% generate input files and store them in datain folder
load('inputFiles');

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
