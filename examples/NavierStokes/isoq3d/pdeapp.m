% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

load('solp2.mat');

% initialize pde structure and mesh structure
[pde,~] = initializeexasim();
pde.model = "ModelD";  
pde.modelfile = "pdemodel";

% Choose computing platform and set number of processors
pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 16;             % number of MPI processors
pde.porder = 2;          % polynomial degree
pde.pgauss = 2*pde.porder;
pde.hybrid = 1;               % 0 -> LDG, 1 -> HDG
pde.debugmode = 0;
pde.nd = 3;

gam = 1.4;                      % specific heat ratio
Re = 1.835e5;                     % Reynolds number
Pr = 0.71;                      % Prandtl number    
Minf = 7;                       % Mach number
Tref  = 124.49;
Twall = 294.44;
pinf = 1/(gam*Minf^2);
Tinf = pinf/(gam-1);
alpha = 0;                % angle of attack
rinf = 1.0;                     % freestream density
ruinf = cos(alpha);             % freestream horizontal velocity
rvinf = sin(alpha);             % freestream vertical velocity
rwinf = 0;
pinf = 1/(gam*Minf^2);          % freestream pressure
rEinf = 0.5+pinf/(gam-1);       % freestream energy

pde.physicsparam = [gam Re Pr Minf rinf ruinf rvinf rwinf rEinf Tinf Tref Twall];
pde.tau = 4.0;                  % DG stabilization parameter
pde.GMRESrestart = 500;         %try 50
pde.GMRESortho = 1;
pde.linearsolvertol = 1e-6; % GMRES tolerance
pde.linearsolveriter = 1000; %try 100
pde.preconditioner = 1;
pde.RBdim = 0;
pde.ppdegree = 0;
pde.NLtol = 1e-6;              % Newton tolerance
pde.NLiter = 10;                % Newton iterations
pde.matvectol=1e-6;             % tolerance for matrix-vector multiplication

pde.dt = [1e-3*(1.15.^(0:11)) 0.005*ones(1,100)];
pde.nstage = 1;
pde.torder = 1;
pde.saveSolFreq = 4;

mesh.dgnodes(:,2,:) = mesh.dgnodes(:,2,:) - 1e-4;
mesh.p(2,:) = mesh.p(2,:) - 1e-4;

mesh3d = mkmesh_isoq3d(pde.porder, 6);
% symmetry, symmetry, outflow, wall, inflow
mesh3d.boundarycondition = [4, 4, 2, 3, 1]; % Set boundary condition for each boundary

xdg = mesh3d.dgnodes(:,1,:);
xdg(:,2,:) = sqrt(mesh3d.dgnodes(:,2,:).^2 + mesh3d.dgnodes(:,3,:).^2 + 1e-12);
costhe = mesh3d.dgnodes(:,2,:)./xdg(:,2,:);
sinthe = mesh3d.dgnodes(:,3,:)./xdg(:,2,:);

udg3d = fieldatuniquedgnodes(mesh, master, sol, xdg);
[npe, nd, ne] = size(mesh3d.dgnodes);
mesh3d.udg = zeros(npe, 5, ne);
mesh3d.udg(:,1:2,:) = udg3d(:,1:2,:); % rho, rho u
mesh3d.udg(:,3,:) = costhe.*udg3d(:,3,:); % rho v
mesh3d.udg(:,4,:) = sinthe.*udg3d(:,3,:); % rho w
mesh3d.udg(:,5,:) = udg3d(:,4,:);
% mesh3d.udg(:,1:2,:) = 1;
% mesh3d.udg(:,5,:) = pde.physicsparam(9);
mesh3d.vdg = fieldatuniquedgnodes(mesh, master, 3*mesh.vdg, xdg);


[pde,mesh,master,dmd] = preprocessing(pde,mesh3d);
% kkgencode(pde);
% compilerstr = cmakecompile(pde); % use cmake to compile C++ source codes 
% runcode(pde, 1); % run C++ code
% sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');




% [mesh2d, xc, xd] = mkmesh_isoq2d(pde.porder);
% mesh2d.vdg = fieldatdgnodes(mesh, master, mesh.vdg, mesh2d.dgnodes); 
% mesh2d.udg = fieldatdgnodes(mesh, master, sol, mesh2d.dgnodes); 
% figure(1); clf; scaplot(mesh2d, eulereval(mesh2d.udg, 'M',gam,Minf),[],1); colorbar; axis equal; axis tight;set(gca,'FontSize',16);
% figure(2); clf; scaplot(mesh2d, mesh2d.vdg(:,1,:),[],1); colorbar; axis equal; axis tight;set(gca,'FontSize',16);
% figure(3); clf; meshplot(mesh,1); hold on;
% x = mesh2d.dgnodes(:,1,:);
% y = mesh2d.dgnodes(:,2,:);
% plot(x(:), y(:), 'o');

% figure(4); clf; meshplot(mesh,1); hold on;
% x = xdg(:,1,:);
% y = xdg(:,2,:);
% plot(x(:), y(:), 'o');
% axis on;

% mesh3d.udg(:,6,:) = mesh3d.vdg;
% master = Master(pde);
% mesh3d.xpe = master.xpe;
% pde.visscalars = {"density", 1, "energy", 5, "av", 6};                % list of scalar fields for visualization
% pde.visvectors = {"velocity", [2 3 4]}; % list of vector fields for visualization
% vis(mesh3d.udg,pde,mesh3d);                        % visualize the numerical solution
% disp("Done!");
