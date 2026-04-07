% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
ii = ii(end);
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde structure and mesh structure
[pde,~] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
pde.platform = "cpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 4;              % number of MPI processors
pde.hybrid = 1;
pde.porder = 2;          % polynomial degree
pde.gencode = 1;

mesh = mkmesh_cylns(pde.porder);
% iso-thermal wall, supersonic outflow, supersonic inflow
mesh.boundarycondition = [3;2;1]; 

gam = 1.4;                       % specific heat ratio
Re = 2.36e5;                     % Reynolds number
Pr = 0.71;                       % Prandtl number    
Minf = 8.98;                     % Mach number
Tref  = 901;
Twall = 300;
pinf = 1/(gam*Minf^2);
Tinf = pinf/(gam-1);
alpha = 0;                % angle of attack
rinf = 1.0;                     % freestream density
ruinf = cos(alpha);             % freestream horizontal velocity
rvinf = sin(alpha);             % freestream vertical velocity
pinf = 1/(gam*Minf^2);          % freestream pressure
rEinf = 0.5+pinf/(gam-1);       % freestream energy

pde.physicsparam = [gam Re Pr Minf rinf ruinf rvinf rEinf Tinf Tref Twall];
pde.tau = 10.0;                  % DG stabilization parameter
pde.GMRESrestart = 500;         %try 50
pde.linearsolvertol = 1e-6; % GMRES tolerance
pde.linearsolveriter = 500; %try 100
pde.RBdim = 0;
pde.ppdegree = 0;
pde.NLtol = 1e-6;              % Newton tolerance
pde.NLiter = 20;                 % Newton iterations

% initial artificial viscosity
mesh.f = facenumbering(mesh.p,mesh.t,pde.elemtype,mesh.boundaryexpr,mesh.periodicexpr);
dist = meshdist3(mesh.f,mesh.dgnodes,mesh.perm,[1]); % distance to the wall
mesh.vdg = zeros(size(mesh.dgnodes,1),2,size(mesh.dgnodes,3));
nm = 30;
mesh.vdg(:,1,:) = 0.06.*tanh(dist*nm);

% intial solution
ui = [rinf ruinf rvinf rEinf];
UDG = initu(mesh,{ui(1),ui(2),ui(3),ui(4),0,0,0,0,0,0,0,0}); % freestream 
UDG(:,2,:) = UDG(:,2,:).*tanh(nm*dist);
UDG(:,3,:) = UDG(:,3,:).*tanh(nm*dist);
TnearWall = Tinf * (Twall/Tref-1) * exp(-nm*dist) + Tinf;
UDG(:,4,:) = TnearWall + 0.5*(UDG(:,2,:).*UDG(:,2,:) + UDG(:,3,:).*UDG(:,3,:));
mesh.udg = UDG;

[sol,pde,mesh,master] = exasim(pde,mesh);
figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],2,1);

disp("Iter 2")
mesh.vdg(:,1,:) = 0.03.*tanh(dist*nm);
mesh.udg = sol;
[pde,mesh,master,dmd] = preprocessing(pde,mesh);
runcode(pde, 1); % run C++ code
sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],2,1);

disp("Iter 3")
mesh.vdg(:,1,:) = 0.02.*tanh(dist*nm);
mesh.udg = sol;
[pde,mesh,master,dmd] = preprocessing(pde,mesh);
runcode(pde, 1); % run C++ code
sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[]);
figure(2); clf; scaplot(mesh, (Tref/Tinf)*eulereval(sol, 't',gam,Minf),[]);

disp("Iter 4")
mesh.vdg(:,1,:) = 0.015.*tanh(dist*nm);
mesh.udg = sol;
[pde,mesh,master,dmd] = preprocessing(pde,mesh);
runcode(pde, 1); % run C++ code
sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],1);

% disp("Iter 5")
% nm = 15; mesh.vdg(:,1,:) = 0.015.*tanh(dist*nm);
% mesh.udg = sol;
% [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% runcode(pde, 1); % run C++ code
% sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],1);

% disp("Iter 6")
% nm = 5; mesh.vdg(:,1,:) = 0.015.*tanh(dist*nm);
% mesh.udg = sol;
% [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% runcode(pde, 1); % run C++ code
% sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],1);
% colorbar; colormap('jet'); axis on; axis equal; axis tight; set(gca,'FontSize',16);


% disp("Iter 5")
% mesh.vdg(:,1,:) = 0.01.*tanh(dist*nm);
% mesh.udg = sol;
% [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% runcode(pde, 1); % run C++ code
% sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],1);
% colorbar; colormap('jet'); axis on; axis equal; axis tight; set(gca,'FontSize',16);
% 
% disp("Iter 6")
% nm = 20;
% mesh.vdg(:,1,:) = 0.01.*tanh(dist*nm);
% mesh.udg = sol;
% [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% runcode(pde, 1); % run C++ code
% sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],1);
% colorbar; colormap('jet'); axis on; axis equal; axis tight; set(gca,'FontSize',16);
% 
% disp("Iter 7")
% nm = 10;
% mesh.vdg(:,1,:) = 0.01.*tanh(dist*nm);
% mesh.udg = sol;
% [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% runcode(pde, 1); % run C++ code
% sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],1);
% colorbar; colormap('jet'); axis on; axis equal; axis tight; set(gca,'FontSize',16);
% 
% disp("Iter 8")
% nm = 5;
% mesh.vdg(:,1,:) = 0.01.*tanh(dist*nm);
% mesh.udg = sol;
% [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% runcode(pde, 1); % run C++ code
% sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],1);
% colorbar; colormap('jet'); axis on; axis equal; axis tight; set(gca,'FontSize',16);
% 
% disp("Iter 9")
% divmax = 5; pdeapp_hm;
% mesh.vdg(:,1,:) = 0.05*av;
% mesh.udg = sol;
% pde.GMRESrestart = 500;         %try 50
% pde.linearsolveriter = 500; %try 100
% preprocessing(pde,mesh);
% runcode(pde, 1); % run C++ code
% sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[0 Minf],1); colorbar;
% figure(2); clf; scaplot(mesh, mesh.vdg(:,1,:),[],2); colorbar;
% 
% disp("Iter 10")
% divmax = 5; pdeapp_hm;
% mesh.vdg(:,1,:) = 0.03*av;
% mesh.udg = sol;
% pde.GMRESrestart = 500;         %try 50
% pde.linearsolveriter = 500; %try 100
% preprocessing(pde,mesh);
% runcode(pde, 1); % run C++ code
% sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[0 Minf],1); colorbar;
% figure(2); clf; scaplot(mesh, mesh.vdg(:,1,:),[],2); colorbar;
% 
% disp("Iter 10")
% divmax = 5; pdeapp_hm;
% mesh.vdg(:,1,:) = 0.02*av;
% mesh.udg = sol;
% pde.GMRESrestart = 500;         %try 50
% pde.linearsolveriter = 500; %try 100
% preprocessing(pde,mesh);
% runcode(pde, 1); % run C++ code
% sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[0 Minf],1); colorbar;
% figure(2); clf; scaplot(mesh, mesh.vdg(:,1,:),[],2); colorbar;
% 
% disp("Iter 11")
% divmax = 5; pdeapp_hm;
% mesh.vdg(:,1,:) = 0.015*av;
% mesh.udg = sol;
% pde.GMRESrestart = 500;         %try 50
% pde.linearsolveriter = 500; %try 100
% preprocessing(pde,mesh);
% runcode(pde, 1); % run C++ code
% sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[0 Minf],1); colorbar;
% figure(2); clf; scaplot(mesh, eulereval(sol, 't',gam,Minf),[],1); colorbar;
% figure(3); clf; scaplot(mesh, mesh.vdg(:,1,:),[],2); colorbar;
% 
% disp("Iter 7")
% divmax = 5; pdeapp_hm;
% mesh.vdg(:,1,:) = 0.024*av;
% mesh.udg = sol;
% preprocessing(pde,mesh);
% runcode(pde, 1); % run C++ code
% sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[0 Minf],1); colorbar;
% figure(2); clf; scaplot(mesh, mesh.vdg(:,1,:),[],2); colorbar;
