% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde structure and mesh structure
pdet = initializepde("v1");

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pdet.model = "ModelC";          % ModelC, ModelD, ModelW
pdet.modelfile = "pdemodel_transport";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
%pdet.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pdet.mpiprocs = 1;              % number of MPI processors
pdet.hybrid = 1;

% Set discretization parameters, physical parameters, and solver parameters
pdet.porder = 3;          % polynomial degree
pdet.pgauss = 2*pdet.porder;
pdet.torder = 2;          % time-stepping order of accuracy
pdet.nstage = 2;          % time-stepping number of stages
pdet.tau = 1.0;               % DG stabilization parameter
pdet.dt = 0.02*ones(1,50);   % time step sizes
pdet.soltime = 1:length(pdet.dt); % steps at which solution are collected
pdet.visdt = 0.02; % visualization timestep size
pdet.physicsparam = [1 theta a1 a2 a];  

% % call exasim to generate and run C++ code to solve the PDE model
[solt,pdet] = exasim(pdet,mesh);

% for i = 1:length(pdet.dt)
%     figure(1); clf; scaplot(mesh,solt(:,1,:,i),[],2); axis on; axis equal; axis tight;
%     pause(0.5);
% end


% mesh.porder = pdet.porder;
% mesh.dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pdet.porder);
% for i = 1:size(sol,4)
%   figure(1); clf; scaplot(mesh,sol(:,1,:,i),[-1 1],2,1); axis on; axis equal; axis tight;
% end
% 
