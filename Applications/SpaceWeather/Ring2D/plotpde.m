% Specify an Exasim version to run
version = "Version0.1";

% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim(version);

% Choose computing platform and set number of processors
% pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors

[pde,mesh] = pdeparams(pde,mesh);

% search compilers and set options
pde = setcompilers(pde);       

% generate input files and store them in datain folder
[pde,mesh,master,dmd] = preprocessing(pde,mesh);

% get solution from output files in dataout folder
sol = fetchsolution(pde,master,dmd);


gam = pde.physicsparam(1);
R = (gam-1)*pde.physicsparam(7)/gam;

nt = length(pde.soltime);
sol2 = zeros(master.npe,5,pde.ne,nt);

for iT = 1:nt
   sol2(:,1:4,:,iT) = sol(:,1:4,:,iT); 
   p = eulereval(sol(:,:,:,iT),'p',pde.physicsparam(1),pde.physicsparam(4));
   rho = eulereval(sol(:,:,:,iT),'r',pde.physicsparam(1),pde.physicsparam(4));
   sol2(:,5,:,iT) = p./(R*rho);
end

% visualize the numerical solution of the PDE model using Paraview
pde.visscalars = {"density", 1, "energy", 4, "temperature", 5};  % list of scalar fields for visualization
pde.visvectors = {"momentum", [2, 3]}; % list of vector fields for visualization
xdg = vis(sol2,pde,mesh); % visualize the numerical solution
disp("Done!");