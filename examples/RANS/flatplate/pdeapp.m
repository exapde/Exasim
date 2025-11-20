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
pde.linearsolveriter = 1000; %try 100
pde.preconditioner = 1;
pde.RBdim = 0;
pde.ppdegree = 0;
pde.NLtol = 1e-8;              % Newton tolerance
pde.NLiter = 3;                % Newton iterations
pde.gencode = 0;

ntime = 11;
dt=1e-2*2.^(0:ntime);
dt=repmat(dt,[2 1]);
pde.dt=dt(:);

% % BJ, ASM, ASM-PP(5), ASM-PP(10)
% ntime = 200;
% % [1e-3 2e-3 4e-3 8e-3 16e-3 32e-3]
% dt = 16e-3*ones(ntime,1);
% pde.dt=dt(:);

mesh = mkmesh_flatplate(1.5,1,1,pde.porder);

master = Master(pde);
dist = meshdist3(mesh.f,mesh.dgnodes,master.perm,3); % distance to the wall

mesh.boundarycondition = [1 1 2 3 1]; 
mesh.vdg = 1e-5*exp(-10000*(mesh.dgnodes(:,1,:).^2 + mesh.dgnodes(:,2,:).^2));  
mesh.vdg(:,2,:) = dist;

% intial solution
Tinf = pinf/(gam-1);
nm = 100;
ui = [rinf ruinf rvinf rEinf rNinf];
UDG = initu(mesh,{ui(1),ui(2),ui(3),ui(4),ui(5),0,0,0,0,0,0,0,0,0,0}); % freestream 
UDG(:,2,:) = UDG(:,2,:).*tanh(nm*dist);
UDG(:,3,:) = UDG(:,3,:).*tanh(nm*dist);
UDG(:,4,:) = Tinf + 0.5*(UDG(:,2,:).*UDG(:,2,:) + UDG(:,3,:).*UDG(:,3,:));
mesh.udg = UDG;

% figure(3); clf; scaplot(mesh,ruinf*tanh(100*dist),[],2); 
% axis off; axis tight; axis square;

[sol,pde,mesh,master,dmd] = exasim(pde,mesh);

% figure(1); clf; scaplot(mesh, eulereval(sol(:,:,:,end), 'M',gam,Minf),[],2);
% figure(2); clf; scaplot(mesh, eulereval(sol(:,:,:,end), 'u',gam,Minf),[],2);
% figure(3); clf; scaplot(mesh, sol(:,5,:,end),[],2);

figure(1); clf; scaplot(mesh, eulereval(sol(:,:,:,end), 'M',gam,Minf),[],2); 
colormap("jet"); axis on; axis tight; axis normal; set(gca,'fontsize', 18);
colorbar("Location", "NorthOutside");

exportgraphics(gca,"ransnacap3mach.png",'Resolution',200);

figure(3); clf; scaplot(mesh, eulereval(sol(:,:,:,end), 'u',gam,Minf),[0 1],2); 
colormap("jet"); axis on; axis tight; axis normal; set(gca,'fontsize', 18);
colorbar("Location", "NorthOutside");
exportgraphics(gca,"ransfpp3uvel.png",'Resolution',200);
 
figure(3); clf; scaplot(mesh, sol(:,5,:,end)*Re,[0 1.8e-4*Re],2);
colormap("jet"); axis on; axis tight; axis normal; set(gca,'fontsize', 18);
colorbar("Location", "NorthOutside");
exportgraphics(gca,"ransfpp3visc.png",'Resolution',200);

