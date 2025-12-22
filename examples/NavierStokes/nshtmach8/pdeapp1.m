% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
ii = ii(end);
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde structure and mesh structure
[pde{1},~] = initializeexasim();
pde{1}.buildpath = string(pwd);

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde{1}.model = "ModelD";          % ModelC, ModelD, ModelW
pde{1}.modelfile = "pdemodel1";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
pde{1}.platform = "cpu";         % choose this option if NVIDIA GPUs are available
pde{1}.mpiprocs = 2;              % number of MPI processors
pde{1}.hybrid = 1;
pde{1}.porder = 4;          % polynomial degree
pde{1}.gencode = 1;
pde{1}.Cxxpreprocessing = 0;

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

pde{1}.physicsparam = [gam Re Pr Minf rinf ruinf rvinf rEinf Tinf Tref Twall];
pde{1}.tau = 1.0;                  % DG stabilization parameter
pde{1}.GMRESrestart = 100;         %try 50
pde{1}.linearsolvertol = 1e-8; % GMRES tolerance
pde{1}.linearsolveriter = 100; %try 100
pde{1}.RBdim = 0;
pde{1}.ppdegree = 20;
pde{1}.NLtol = 1e-6;              % Newton tolerance
pde{1}.NLiter = 30;                 % Newton iterations
pde{1}.matvectol=1e-6;             % tolerance for matrix-vector multiplication

mesh{1} = mkmesh_cylns(pde{1}.porder);
% iso-thermal wall, supersonic outflow, supersonic inflow
mesh{1}.boundarycondition = [3;2;1]; 

% initial artificial viscosity
mesh{1}.f = facenumbering(mesh{1}.p,mesh{1}.t,pde{1}.elemtype,mesh{1}.boundaryexpr,mesh{1}.periodicexpr);
dist = meshdist3(mesh{1}.f,mesh{1}.dgnodes,mesh{1}.perm,[1]); % distance to the wall
mesh{1}.vdg = zeros(size(mesh{1}.dgnodes,1),2,size(mesh{1}.dgnodes,3));
mesh{1}.vdg(:,1,:) = 0.06.*tanh(dist*30);

% intial solution
ui = [rinf ruinf rvinf rEinf];
UDG = initu(mesh{1},{ui(1),ui(2),ui(3),ui(4),0,0,0,0,0,0,0,0}); % freestream 
UDG(:,2,:) = UDG(:,2,:).*tanh(10*dist);
UDG(:,3,:) = UDG(:,3,:).*tanh(10*dist);
TnearWall = Tinf * (Twall/Tref-1) * exp(-10*dist) + Tinf;
UDG(:,4,:) = TnearWall + 0.5*(UDG(:,2,:).*UDG(:,2,:) + UDG(:,3,:).*UDG(:,3,:));
mesh{1}.udg = UDG;

% mesh{1}.vdg(:,1,:) = 0.025.*tanh(dist*5);
% mesh{1}.udg = solns;

[sol{1},pde{1},mesh{1}] = exasim(pde{1},mesh{1});

disp("Iter 2")
mesh{1}.vdg(:,1,:) = 0.04.*tanh(dist*30);
mesh{1}.udg = sol{1};
[pde{1},mesh{1},master{1},dmd{1}] = preprocessing(pde{1},mesh{1});
runcode(pde{1}, 1); % run C++ code
sol{1} = fetchsolution(pde{1},master{1},dmd{1}, pde{1}.buildpath + '/dataout');

disp("Iter 3")
mesh{1}.vdg(:,1,:) = 0.025.*tanh(dist*30);
mesh{1}.udg = sol{1};
[pde{1},mesh{1},master{1},dmd{1}] = preprocessing(pde{1},mesh{1});
runcode(pde{1}, 1); % run C++ code
sol{1} = fetchsolution(pde{1},master{1},dmd{1}, pde{1}.buildpath + '/dataout');

disp("Iter 4")
mesh{1}.vdg(:,1,:) = 0.025.*tanh(dist*5);
mesh{1}.udg = sol{1};
[pde{1},mesh{1},master{1},dmd{1}] = preprocessing(pde{1},mesh{1});
runcode(pde{1}, 1); % run C++ code
sol{1} = fetchsolution(pde{1},master{1},dmd{1}, pde{1}.buildpath + '/dataout');
mesh{1}.udg = sol{1};

figure(1); clf; scaplot(mesh{1}, eulereval(sol{1}, 'M',gam,Minf),[]);
figure(2); clf; scaplot(mesh{1}, (Tref/Tinf)*eulereval(sol{1}, 't',gam,Minf),[]);




