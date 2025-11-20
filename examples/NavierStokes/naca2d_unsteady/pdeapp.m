% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

porder = 3;                     % polynomial degree
gam = 1.4;                      % gas constant
Minf = 0.25;                    % freestream mach number
tau = 1;                        % stabilization parameter
alpha = 8*pi/180;               % angle of attack
rinf = 1.0;                     % freestream density
ruinf = cos(alpha);             % freestream horizontal velocity
rvinf = sin(alpha);             % freestream vertical velocity
pinf = 1/(gam*Minf^2);          % freestream pressure
rEinf = 0.5+pinf/(gam-1);       % freestream energy
Re = 5000;                       % Reynolds number 
Pr = 0.72;                      % Prandtl number 
ui = [ 1, cos(alpha), sin(alpha), 0.5+pinf/(gam-1)];

% initialize pde structure and mesh structure
[pde,~] = initializeexasim();

pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
pde.platform = "cpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 8;              % number of MPI processors
pde.hybrid = 1;
pde.debugmode = 0;
pde.porder = porder;
pde.pgauss = 2*porder;

pde.physicsparam = [gam Re Pr Minf rinf ruinf rvinf rEinf];
pde.tau = tau;              % DG stabilization parameter
pde.GMRESrestart = 100;
pde.GMRESortho = 1;
pde.linearsolvertol = 1e-6; % GMRES tolerance
pde.linearsolveriter = 100;
pde.preconditioner = 1;
pde.NLtol = 1e-8;
pde.ppdegree = 1;
pde.RBdim = 0;
pde.gencode = 1;

pde.torder = 2;
pde.nstage = 2;
pde.dt = [0.005 0.01 0.02 0.04 0.06 0.08*ones(1,195)];
pde.saveSolFreq = 4;

% naca mesh
mesh = mkmesh_naca0012(porder,1,3);

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh] = exasim(pde,mesh);


% % plot solution
mesh.porder = porder;
mesh.dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,porder);    

for i = 1:size(sol,4)
    figure(1); clf; scaplot(mesh,eulereval(sol(:,1:4,:,i),'u',gam),[],1);
    colorbar; axis tight; axis equal; axis([-0.5 4 -0.5 1]);
    pause(1);
end

