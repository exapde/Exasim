% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
ii = ii(end);
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel_avcont";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors
pde.hybrid = 1;
% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;          % polynomial degree
% pde.torder = 1;          % time-stepping order of accuracy
% pde.nstage = 1;          % time-stepping number of stages
% pde.dt = 1e-4*ones(1,10000)*5;   % time step sizes
% pde.visdt = pde.dt(1);         % visualization timestep size
% pde.saveSolFreq = 100;          % solution is saved every 10 time steps
% pde.soltime = 100:100:length(pde.dt); % steps at which solution are collected

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
pde.GMRESrestart=39;            % number of GMRES restarts
pde.linearsolvertol=1e-2;     % GMRES tolerance
pde.linearsolveriter=80;        % number of GMRES iterations
pde.GMRESrestart = 400;
pde.linearsolvertol = 1e-3; % GMRES tolerance
pde.linearsolveriter = 400;

pde.ppdegree = 20;
pde.precMatrixType=2;           % preconditioning type
pde.ptcMatrixType=0;
pde.NLtol = 1e-12;              % Newton tolerance
pde.NLiter = 30;                 % Newton iterations
pde.matvectol=1e-6;             % tolerance for matrix-vector multiplication
pde.timestepOffset=0;
pde.gencode = 1;
pde.AV = 1;

% [mesh.p,mesh.t,mesh.dgnodes] = mkmesh_circincirc_Ma17b(pde.porder,201,201,1,3,4);
mesh = mkmesh_square(51,21,pde.porder,1,1,1,1,1);
mesh.p(1,:) = logdec(mesh.p(1,:), 4.5);
mesh.dgnodes(:,1,:) = logdec(mesh.dgnodes(:,1,:), 4.5);
mesh = mkmesh_halfcircle(mesh, 1, 3, 4.7, pi/2, 3*pi/2);
mesh.porder = pde.porder;
mesh.boundaryexpr = {@(p) sqrt(p(1,:).^2+p(2,:).^2)<1+1e-6, @(p) p(1,:)>-1e-7, @(p) abs(p(1,:))<20};
mesh.periodicexpr = {};
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
mesh.vdg(:,1,:) = 0.06.*tanh(dist*30);
% mesh.dgnodes(:,3,:) = 0.06.*tanh(dist*30);

% intial solution
ui = [rinf ruinf rvinf rEinf];
UDG = initu(mesh,{ui(1),ui(2),ui(3),ui(4),0,0,0,0,0,0,0,0}); % freestream 

UDG(:,2,:) = UDG(:,2,:).*tanh(10*dist);
UDG(:,3,:) = UDG(:,3,:).*tanh(10*dist);
TnearWall = Tinf * (Twall/Tref-1) * exp(-10*dist) + Tinf;
UDG(:,4,:) = TnearWall + 0.5*(UDG(:,2,:).*UDG(:,2,:) + UDG(:,3,:).*UDG(:,3,:));
mesh.udg = UDG;


% search compilers and set options
pde = setcompilers(pde);       

% generate input files and store them in datain folder
[pde,mesh,master,dmd] = preprocessing(pde,mesh);

% generate source codes and store them in app folder
if pde.gencode==1
  %gencode(pde);
  kkgencode(pde);
end

compilerstr = cmakecompile(pde); % use cmake to compile C++ source codes 

runcode(pde, 1); % run C++ code

% get solution from output files in dataout folder
sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],2,1);

disp("Iter 2")
mesh.vdg(:,1,:) = 0.04.*tanh(dist*30);
mesh.udg = sol;
[pde,mesh,master,dmd] = preprocessing(pde,mesh);
runcode(pde, 1); % run C++ code
sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],2,1);

disp("Iter 3")
mesh.vdg(:,1,:) = 0.035.*tanh(dist*30);
mesh.udg = sol;
[pde,mesh,master,dmd] = preprocessing(pde,mesh);
runcode(pde, 1); % run C++ code
sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],2,1);
sol = real(sol);

disp("Iter 4")
mesh.vdg(:,1,:) = 0.03.*tanh(dist*30);
mesh.udg = sol;
[pde,mesh,master,dmd] = preprocessing(pde,mesh);
runcode(pde, 1); % run C++ code
sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[]);
UDG0 = sol;

%%
disp("Iter 5")
UDG = UDG0;
disp("Iter 5")
mesh.vdg(:,1,:) = 0.04.*tanh(dist*20);
mesh.udg = UDG;
[pde,mesh,master,dmd] = preprocessing(pde,mesh);
runcode(pde, 1); % run C++ code
sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[]);
UDG = sol;

% mesh.vdg(:,1,:) = 0.04.*tanh(dist*20);
% mesh.udg = UDG;
% [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% runcode(pde, 1); % run C++ code
% sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[]);
% UDG = sol;
% 
% mesh.vdg(:,1,:) = 0.03.*tanh(dist*18);
% mesh.udg = UDG;
% [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% runcode(pde, 1); % run C++ code
% sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[]);
% UDG = sol;
% 
% mesh.vdg(:,1,:) = 0.03.*tanh(dist*16);
% mesh.udg = UDG;
% [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% runcode(pde, 1); % run C++ code
% sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[]);
% UDG = sol;
% 
% mesh.vdg(:,1,:) = 0.03.*tanh(dist*15);
% mesh.udg = UDG;
% [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% runcode(pde, 1); % run C++ code
% sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[]);
% UDG = sol;
% 
% mesh.vdg(:,1,:) = 0.03.*tanh(dist*14);
% mesh.udg = UDG;
% [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% runcode(pde, 1); % run C++ code
% sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[]);
% UDG = sol;
% 
% disp("Iter 5")
% mesh.vdg(:,1,:) = 0.03.*tanh(dist*12);
% mesh.udg = UDG;
% [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% runcode(pde, 1); % run C++ code
% sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[]);
% UDG = sol;
% 
% disp("Iter 5")
% mesh.vdg(:,1,:) = 0.03.*tanh(dist*10);
% mesh.udg = UDG;
% [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% runcode(pde, 1); % run C++ code
% sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[]);
% UDG = sol;
% 
% %%
% disp("Iter 5")
% mesh.vdg(:,1,:) = 0.04.*tanh(dist*8);
% mesh.udg = UDG;
% [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% runcode(pde, 1); % run C++ code
% sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[]);
% UDG = sol;
% 
% disp("Iter 5")
% mesh.vdg(:,1,:) = 0.04.*tanh(dist*5);
% mesh.udg = UDG;
% [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% runcode(pde, 1); % run C++ code
% sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[]);
% UDG = sol;
% 
% 
% dd = [20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5];
% % for i = 1:5
% %     mesh.vdg(:,1,:) = 0.05.*tanh(dist*dd(i));
% %     mesh.udg = UDG;
% %     [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% %     runcode(pde, 1); % run C++ code
% %     sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% %     figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[]);
% %     UDG = sol;
% % end
% % 
% % %%
% % dd = [20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5];
% % for i = 1:20
% %     mesh.vdg(:,1,:) = 0.04.*tanh(dist*dd(i));
% %     mesh.udg = UDG;
% %     [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% %     runcode(pde, 1); % run C++ code
% %     sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% %     figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[]);
% %     UDG = sol;
% % end
% 
% %%
% pde.S0=0.01; pde.lambda = 0.05; pde.kappa=8;
% mesh.dist = dist;
% a = avf(mesh, master, pde, real(UDG_curr));
% figure(1); clf; scaplot(mesh, a)
% mesh.vdg(:,1,:) = a;
% mesh.udg = UDG_curr;
% [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% runcode(pde, 1); % run C++ code
% 
% % 
% % %%
% % % disp("Iter 5")
% % % mesh.vdg(:,1,:) = 0.025.*tanh(dist*30);
% % % mesh.udg = sol;
% % % [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% % % runcode(pde, 1); % run C++ code
% % % sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% % % figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],2,1);
% % % 
% % % disp("Iter 5")
% % % mesh.vdg(:,1,:) = 0.03.*tanh(dist*2.5);
% % % mesh.udg = sol;
% % % [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% % % runcode(pde, 1); % run C++ code
% % % sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% % % figure(1); clf; scaplot(mesh, eulereval(real(sol), 'M',gam,Minf),[],2,1);
% % 
% % % % call exasim to generate and run C++ code to solve the PDE model
% % % [sol,pde,mesh] = exasim(pde,mesh);
% % 
% % % % search compilers and set options
% % % pde = setcompilers(pde);       
% % % 
% % % % generate input files and store them in datain folder
% % % [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% % % 
% % % % generate source codes and store them in app folder
% % % gencode(pde);
% % % 
% % % % compile source codes to build an executable file and store it in app folder
% % % compilerstr = compilecode(pde);
% % % 
% % % % run code
% % % runstr = runcode(pde);
% % % 
% % % % get solution from output files in dataout folder
% % % sol = fetchsolution(pde,master,dmd);
% % % 
% % % % % visualize the numerical solution of the PDE model using Paraview
% % % pde.visscalars = {"density", 1, "energy", 4};  % list of scalar fields for visualization
% % % pde.visvectors = {"momentum", [2, 3]}; % list of vector fields for visualization
% % % xdg = vis(sol,pde,mesh); % visualize the numerical solution
% % % disp("Done!");
