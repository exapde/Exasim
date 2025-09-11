% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde structure and mesh structure
pde = initializeexasim();
pde.model = "ModelD";  
pde.modelfile = "pdemodel";

% Choose computing platform and set number of processors
pde.platform = "cpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;             % number of MPI processors
pde.porder = 3;          % polynomial degree
pde.pgauss = 2*pde.porder;
pde.hybrid = 1;               % 0 -> LDG, 1 -> HDG
pde.debugmode = 0;
pde.nd = 2;

gam = 1.4;
Minf = 0.2;                  % Infinity conditions
Re = 2.1854e+006;
Pr = 0.72;
Tinf = pinf/(gam-1);
alpha = 0;                      % angle of attack
rinf = 1.0;                     % freestream density
ruinf = cos(alpha);             % freestream horizontal velocity
rvinf = sin(alpha);             % freestream vertical velocity
pinf = 1/(gam*Minf^2);          % freestream pressure
rEinf = 0.5+pinf/(gam-1);       % freestream energy
rNinf = 0.25/Re;

pde.physicsparam = [gam Re Pr Minf rinf ruinf rvinf rEinf rNinf];
pde.tau = 2.0;                  % DG stabilization parameter
pde.GMRESrestart = 250;         %try 50
pde.GMRESortho = 1;
pde.linearsolvertol = 1e-6; % GMRES tolerance
pde.linearsolveriter = 500; %try 100
pde.preconditioner = 1;
pde.RBdim = 0;
pde.ppdegree = 0;
pde.NLtol = 1e-7;              % Newton tolerance
pde.NLiter = 4;                % Newton iterations

ntime = 9;
dt=1e-2*2.^(0:ntime);
dt=repmat(dt,[2 1]);
pde.dt=dt(:);

mesh = mkmesh_flatplate(1.5,1,1,pde.porder);

dist = meshdist3(mesh.f,mesh.dgnodes,master.perm,3); % distance to the wall

mesh.boundarycondition = [1 1 2 3 1]; 
mesh.vdg = 0*mesh.dgnodes(:,1,:);  
mesh.vdg(:,2,:) = dist;

[sol,pde,mesh,master,dmd] = exasim(pde,mesh);

figure(1); clf; scaplot(mesh, eulereval(sol(:,:,:,end), 'M',gam,Minf),[],2);
figure(2); clf; scaplot(mesh, eulereval(sol(:,:,:,end), 'u',gam,Minf),[],2);
figure(3); clf; scaplot(mesh, sol(:,5,:,end),[],2);

