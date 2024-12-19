% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
ii = ii(end);
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde structure and mesh structure
[pde,~] = initializeexasim();
pde.buildpath = string(pwd()) + "/ns";

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel_ns";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
pde.platform = "cpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors
pde.hybrid = 1;
pde.porder = 4;          % polynomial degree
pde.gencode = 1;

mesh = mkmesh_cylns(pde.porder);
% iso-thermal wall, supersonic outflow, supersonic inflow
mesh.boundarycondition = [3;2;1]; 

gam = 1.4;                      % specific heat ratio
Re = 1.835e5;                     % Reynolds number
Pr = 0.71;                      % Prandtl number    
Minf = 8.03;                     % Mach number
Tref  = 124.49;
Twall = 294.44;
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
pde.tau = 1.0;                  % DG stabilization parameter
pde.GMRESrestart = 100;         %try 50
pde.linearsolvertol = 1e-8; % GMRES tolerance
pde.linearsolveriter = 100; %try 100
pde.RBdim = 0;
pde.ppdegree = 20;
pde.NLtol = 1e-6;              % Newton tolerance
pde.NLiter = 30;                 % Newton iterations
pde.matvectol=1e-6;             % tolerance for matrix-vector multiplication

% initial artificial viscosity
mesh.f = facenumbering(mesh.p,mesh.t,pde.elemtype,mesh.boundaryexpr,mesh.periodicexpr);
dist = meshdist3(mesh.f,mesh.dgnodes,mesh.perm,[1]); % distance to the wall
mesh.vdg = zeros(size(mesh.dgnodes,1),2,size(mesh.dgnodes,3));
mesh.vdg(:,1,:) = 0.06.*tanh(dist*30);

% intial solution
ui = [rinf ruinf rvinf rEinf];
UDG = initu(mesh,{ui(1),ui(2),ui(3),ui(4),0,0,0,0,0,0,0,0}); % freestream 
UDG(:,2,:) = UDG(:,2,:).*tanh(10*dist);
UDG(:,3,:) = UDG(:,3,:).*tanh(10*dist);
TnearWall = Tinf * (Twall/Tref-1) * exp(-10*dist) + Tinf;
UDG(:,4,:) = TnearWall + 0.5*(UDG(:,2,:).*UDG(:,2,:) + UDG(:,3,:).*UDG(:,3,:));
mesh.udg = UDG;

[sol,pde,mesh,master] = exasim(pde,mesh);
figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],2,1);

disp("Iter 2")
mesh.vdg(:,1,:) = 0.04.*tanh(dist*30);
mesh.udg = sol;
[pde,mesh,master,dmd] = preprocessing(pde,mesh);
runcode(pde, 1); % run C++ code
sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],2,1);

disp("Iter 3")
mesh.vdg(:,1,:) = 0.025.*tanh(dist*30);
mesh.udg = sol;
[pde,mesh,master,dmd] = preprocessing(pde,mesh);
runcode(pde, 1); % run C++ code
sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],2,1);

disp("Iter 4")
mesh.vdg(:,1,:) = 0.025.*tanh(dist*5);
mesh.udg = sol;
[pde,mesh,master,dmd] = preprocessing(pde,mesh);
runcode(pde, 1); % run C++ code
sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[]);
figure(2); clf; scaplot(mesh, (Tref/Tinf)*eulereval(sol, 't',gam,Minf),[]);

% u = getsolutiononboundary(sol, mesh.f, mesh.perm, 1);
% r  = u(:,1,:);
% uv = u(:,2,:)./u(:,1,:);
% vv = u(:,3,:)./u(:,1,:);
% p  = (gam-1)*(u(:,4,:) - 0.5*(u(:,2,:).*uv + u(:,3,:).*vv));
% T = p./((gam-1)*r);

% [hfx, hfy, Tx, Ty, kappa] = heatflux(sol(:,1:4,:), sol(:,5:12,:), pde.physicsparam);
% figure(3); clf; scaplot(mesh, hfx);
% figure(4); clf; scaplot(mesh, hfy);
% 
% figure(5); clf; scaplot(mesh, Tx);
% figure(6); clf; scaplot(mesh, Ty);
% figure(7); clf; scaplot(mesh, kappa);
% 
% fileID = fopen(pde.buildpath + "/dataout/out_uhat_np0.bin",'r');
% UH = fread(fileID,'double');
% fclose(fileID);
% UH = reshape(UH,4,[]);    
% 
% ibwall = 1;
% UDGb = getsolutiononboundary(sol, mesh.f, mesh.perm, ibwall);
% UHb = getuhonboundary(UH, dmd{1}.elemcon, mesh.f, mesh.perm, ibwall);
% THb = eulereval(UHb, 't',gam,Minf);
% TUb = eulereval(UDGb, 't',gam,Minf);
% 
% [hfbx, hfby] = heatflux(UHb, UDGb(:,5:12,:), pde.physicsparam);
% hfbx = squeeze(hfbx);
% hfby = squeeze(hfby);
% hfbd = squeeze(pde.tau*(UDGb(:,4,:) - UHb(:,4,:)));
% 
% %[sqrt(sum(hfbx.*hfbx)); sqrt(sum(hfby.*hfby)); sqrt(sum(hfbd.*hfbd))]
% 

