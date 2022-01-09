% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
%pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 2;          % polynomial degree
pde.torder = 1;          % time-stepping order of accuracy
pde.nstage = 1;          % time-stepping number of stages
pde.dt = 1e-4*ones(1,10000)*5;   % time step sizes
pde.visdt = pde.dt(1);         % visualization timestep size
pde.saveSolFreq = 100;          % solution is saved every 10 time steps
pde.soltime = 100:100:length(pde.dt); % steps at which solution are collected

gam = 1.4;                      % specific heat ratio
Re = 376930;                     % Reynolds number
Pr = 0.71;                      % Prandtl number    
Minf = 17.605;                     % Mach number
Tref  = 200.;
Twall = 500.;
pinf = 1/(gam*Minf^2);
Tinf = pinf/(gam-1);
alpha = 0;                % angle of attack
rinf = 1.0;                     % freestream density
ruinf = cos(alpha);             % freestream horizontal velocity
rvinf = sin(alpha);             % freestream vertical velocity
pinf = 1/(gam*Minf^2);          % freestream pressure
rEinf = 0.5+pinf/(gam-1);       % freestream energy
avb = 10.0;                     % bulk viscosity parameter
avk = 0.5;                      % thermal viscosity parameter 
avs = 0.0;                      % shear viscosity parameter
sb0 = 0.02;                     % cutoff  dilatation
sb1 = 2.5;                      % maximum dilatation 
pde.physicsparam = [gam Re Pr Minf rinf ruinf rvinf rEinf Tinf Tref Twall avb avk avs pde.porder sb0 sb1];
pde.tau = 4.0;                  % DG stabilization parameter
pde.GMRESrestart=29;            % number of GMRES restarts
pde.linearsolvertol=1e-8;     % GMRES tolerance
pde.linearsolveriter=30;        % number of GMRES iterations
pde.precMatrixType=2;           % preconditioning type
pde.ptcMatrixType=0;
pde.NLtol = 1e-11;              % Newton tolerance
pde.NLiter = 3;                 % Newton iterations
pde.matvectol=1e-6;             % tolerance for matrix-vector multiplication
pde.AV = 1;                     % artificial viscosity
pde.timestepOffset=0;

[mesh.p,mesh.t,mesh.dgnodes] = mkmesh_circincirc_Ma17b(pde.porder,201,201,1,3,4);
mesh.boundaryexpr = {@(p) sqrt(p(1,:).^2+p(2,:).^2)<1+1e-6, @(p) p(1,:)>-1e-7, @(p) abs(p(1,:))<20};
% iso-thermal wall, supersonic outflow, supersonic inflow
mesh.boundarycondition = [3;6;5]; 

% mesh size
pde.pgauss = 2*(pde.porder);
pde.nd = 2;
pde.elemtype = 1;
master = Master(pde);
[~, ~, jac] = volgeom(master.shapent,permute(mesh.dgnodes,[1 3 2]));
hsz = reshape(sqrt(jac),[],1,size(mesh.dgnodes,3));
[~,cgelcon,rowent2elem,colent2elem,~] = mkcgent2dgent(mesh.dgnodes,1e-8);
hh = dg2cg2(max(hsz,0e-5), cgelcon, colent2elem, rowent2elem);
hh = dg2cg2(hh, cgelcon, colent2elem, rowent2elem);

% distance to the wall
mesh.f = facenumbering(mesh.p,mesh.t,pde.elemtype,mesh.boundaryexpr,mesh.periodicexpr);
f = mkf(mesh.t,mesh.f,2);
dist = meshdist(f,mesh.dgnodes,master.perm,[1]); % distance to the wall

mesh.vdg = zeros(size(mesh.dgnodes,1),2,size(mesh.dgnodes,3));
mesh.vdg(:,2,:) = hh.*tanh(1000*dist);

% intial solution
ui = [rinf ruinf rvinf rEinf];
UDG = initu(mesh,{ui(1),ui(2),ui(3),ui(4),0,0,0,0,0,0,0,0}); % freestream 

UDG(:,2,:) = UDG(:,2,:).*tanh(10*dist);
UDG(:,3,:) = UDG(:,3,:).*tanh(10*dist);
TnearWall = Tinf * (Twall/Tref-1) * exp(-10*dist) + Tinf;
UDG(:,4,:) = TnearWall + 0.5*(UDG(:,2,:).*UDG(:,2,:) + UDG(:,3,:).*UDG(:,3,:));
mesh.udg = UDG;

% % call exasim to generate and run C++ code to solve the PDE model
% [sol,pde,mesh] = exasim(pde,mesh);

% search compilers and set options
pde = setcompilers(pde);       

% generate input files and store them in datain folder
[pde,mesh,master,dmd] = preprocessing(pde,mesh);

% generate source codes and store them in app folder
gencode(pde);

% compile source codes to build an executable file and store it in app folder
compilerstr = compilecode(pde);

% run code
runstr = runcode(pde);

% get solution from output files in dataout folder
sol = fetchsolution(pde,master,dmd);

% % visualize the numerical solution of the PDE model using Paraview
pde.visscalars = {"density", 1, "energy", 4};  % list of scalar fields for visualization
pde.visvectors = {"momentum", [2, 3]}; % list of vector fields for visualization
xdg = vis(sol,pde,mesh); % visualize the numerical solution
disp("Done!");
