% clear exasim data from memory
clear pde mesh master dmd sol;

% Specify an Exasim version to run
version = "Version0.3";

% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% create pde and mesh for each PDE model
pdeapp1;
pdeapp2;

% call exasim to generate and run C++ code to solve the PDE models
[sol,pde,mesh,master,dmd,compilerstr,runstr] = exasim(pde,mesh);

% visualize the numerical solution of the PDE model using Paraview
% for m = 1:length(pde)
%     pde{m}.visscalars = {"temperature", 1};  % list of scalar fields for visualization
%     pde{m}.visvectors = {"temperature gradient", [2 3]}; % list of vector fields for visualization
%     pde{m}.visfilename = "dataout" + num2str(m) + "/output";  
%     vis(sol{m},pde{m},mesh{m}); % visualize the numerical solution
% end
% disp("Done!");


