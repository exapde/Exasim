% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
%pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors
pde.hybrid = 1;        

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;          % polynomial degree
pde.physicsparam = [1 100 100];    % unit thermal conductivity
pde.tau = 1.0;           % DG stabilization parameter
pde.GMRESrestart = 50;
% pde.ppdegree = 20;
pde.linearsolvertol = 1e-6;

% create a grid of 8 by 8 on the unit square
[mesh.p,mesh.t] = squaremesh(64,64,1,1);
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(2,:))<1e-8, @(p) abs(p(1,:)-1)<1e-8, @(p) abs(p(2,:)-1)<1e-8, @(p) abs(p(1,:))<1e-8};
mesh.boundarycondition = [1;1;1;1]; % Set boundary condition for each boundary

pde.platform = "gpu";
pde.GMRESortho = 1;
pde.ppdegree = 20;
[sol,pde,mesh] = exasim(pde,mesh);
return;

platforms  = ["cpu", "gpu"];
GMRESortho = [0,1];
ppdegree = 0:20;

for i = 1:2
   for j = 1:2
       for k = 1:21
           pde.platform = platforms(i);         % choose this option if NVIDIA GPUs are available
           pde.GMRESortho = GMRESortho(j);
           pde.ppdegree = ppdegree(k);


           % Create a file name using sprintf
           filename = sprintf('output/platform_%s_GMRESortho_%d_ppdegree_%d.txt', pde.platform, pde.GMRESortho, pde.ppdegree);
           disp(filename)


           diary(filename)
           diary on


           % call exasim to generate and run C++ code to solve the PDE model
           [sol,pde,mesh] = exasim(pde,mesh);


       end
   end
end


% % call exasim to generate and run C++ code to solve the PDE model
% [sol,pde,mesh] = exasim(pde,mesh);

% % visualize the numerical solution of the PDE model using Paraview
% pde.visscalars = {"temperature", 1};  % list of scalar fields for visualization
% pde.visvectors = {"temperature gradient", [2 3]}; % list of vector fields for visualization
% xdg = vis(sol,pde,mesh); % visualize the numerical solution
% disp("Done!");

% mesh.porder = pde.porder;
% mesh.dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);
% figure(1); clf; scaplot(mesh,sol(:,1,:),[-1 1],2,1); axis on; axis equal; axis tight;

        
