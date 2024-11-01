% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
ii = ii(end);
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors
pde.hybrid = 1;
% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 4;          % polynomial degree
porder = pde.porder;
ns=5;
pde.tau = 5.0;                  % DG stabilization parameter
pde.read_uh = 0;
rho_inf = [1.0055500000000001e-09 0.00035318161606 1.5872374699999997e-05 0.00116691299088 1.103201281e-05];
rhou_inf = 9.213932;
rhov_inf = 0.0;
rhoE_inf = 33563.20282790763; %TODO...do i need to change these to be consistent? 
% pde.debugmode=1;

% Reference quantities for nondimensionalizing
rho_ref = 0.001547;  % sum of inflow densities
u_inf = 5976;        % inflow speed
% u_ref = 6.655369622053456e+02;     % inflow speed of sound; used for nondim
u_ref = u_inf;
rhoE_ref = rho_ref * u_ref^2;
T_ref = 901; 
gamma_ref = 1.415845062661413;
mu_ref =    3.8196e-05;
kappa_ref = 0.076685370133649;
% Nondimensional quantities 
cp = 1182.1920097928833;
L_ref = 0.045;
Re = rho_ref * L_ref * u_ref / mu_ref; % Should this be 1? 
Pr = mu_ref * cp / kappa_ref;
Ec = u_ref^2 / (cp * T_ref*(gamma_ref-1));

% Load into Exasim data structures
pde.physicsparam(1:ns) = rho_inf;
pde.physicsparam(ns+1) = rhou_inf;
pde.physicsparam(ns+2) = rhov_inf;
pde.physicsparam(ns+3) = rhoE_inf;
pde.physicsparam(15) = pde.porder;
% pde.physicsparam(16) = 0;
% pde.physicsparam(17) = sb0;
% pde.physicsparam(18) = sbmax;
pde.physicsparam(19) = 1.0e13; %alpha
pde.physicsparam(20) = 1.0e4; %sigma
pde.physicsparam(21) = Pr;
pde.physicsparam(22) = Re;
pde.physicsparam(23) = Ec*0 +1; %avk


pde.externalparam = zeros(1,8);
pde.externalparam(1) = rho_ref;  % rho_inf
pde.externalparam(2) = u_ref;    % u_inf
pde.externalparam(3) = rhoE_ref;  % rhoE_inf
pde.externalparam(4) = T_ref;
pde.externalparam(5) = mu_ref;
pde.externalparam(6) = kappa_ref;
pde.externalparam(7) = cp; % cp_inf
pde.externalparam(8) = L_ref;

pde.GMRESrestart = 50;         %try 50
pde.linearsolvertol = 1e-8; % GMRES tolerance
pde.linearsolveriter = 100; %try 100
pde.RBdim = 0;
pde.ppdegree = 0;
pde.precMatrixType=2;           % preconditioning type
pde.ptcMatrixType=0;
pde.NLtol = 1e-6;              % Newton tolerance
pde.NLiter = 30;                 % Newton iterations
pde.matvectol=1e-6;             % tolerance for matrix-vector multiplication
pde.timestepOffset=0;
pde.dae_alpha = 0.0;
pde.dae_beta = 0.0;
pde.AV = 1;%

% [mesh.p,mesh.t,mesh.dgnodes] = mkmesh_circincirc_Ma17b(pde.porder,201,201,1,3,4);
mesh = mkmesh_square(61,61,pde.porder,1,1,1,1,1);
mesh.p(1,:) = logdec(mesh.p(1,:),0.5);
mesh.dgnodes(:,1,:) = logdec(mesh.dgnodes(:,1,:),0.5);
mesh = mkmesh_halfcircle(mesh,1,2,3,pi/2,3*pi/2);

mesh.porder = pde.porder;
mesh.boundaryexpr = {@(p) sqrt(p(1,:).^2+p(2,:).^2)<1+1e-6, @(p) p(1,:)>-1e-7, @(p) abs(p(1,:))<20};
mesh.periodicexpr = {};
% iso-thermal wall, supersonic outflow, supersonic inflow
mesh.boundarycondition = [3;6;5]; 
pde.AVsmoothingIter = 0;
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

% intial solution
ui = [ rho_inf / rho_ref, rhou_inf / (rho_ref *u_ref), 0.0, rhoE_inf / (rhoE_ref)];
UDG = initu(mesh,{ui(1),ui(2),ui(3),ui(4),ui(5),ui(6),ui(7),ui(8),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0});
UDG(:,ns+1,:) = UDG(:,ns+1,:).*tanh(dist*4);
UDG(:,ns+2,:) = UDG(:,ns+2,:).*tanh(dist*4);
% TnearWall = pinf/(gam-1); % Tinf * (Twall/Tref-1) * exp(-10*dist) + Tinf;
rhoe_dim = rhoE_inf - 1/2 * sum(rho_inf) * u_inf^2;
rhoe = rhoe_dim / (rhoE_ref);
UDG(:,ns+3,:) = rhoe + 1/2 * (UDG(:,ns+1,:).^2 + UDG(:,ns+2,:).^2) ./ sum( UDG(:,1:ns,:), 2 ); 
mesh.udg = UDG;
rhoe_dim = rhoE_inf - 1/2 * sum(rho_inf) * u_inf^2;
rhoe = rhoe_dim / (rhoE_ref);
% mesh.udg = UDG;
WDG = ones(size(mesh.dgnodes,1),1,size(mesh.dgnodes,3));
%%
for i = 1:1
 WDG(:,i,:) = 901.0 / pde.externalparam(4);
end

mesh.udg = UDG;
% WDG = Tsmooth / pde.externalparam(4);
mesh.wdg = WDG;%
mesh.vdg = zeros(size(mesh.dgnodes,1),2,size(mesh.dgnodes,3));
mesh.vdg(:,2,:) = hh.*tanh(1000*dist);
mesh.vdg(:,1,:) = 0.1*tanh(dist*25);
% mesh.dgnodes(:,3,:) = 0.06.*tanh(dist*30);

pde.dae_alpha = 0.0;
pde.dae_beta = 0.0;


% search compilers and set options
pde = setcompilers(pde);       

% generate input files and store them in datain folder
[pde,mesh,master,dmd] = preprocessing(pde,mesh);

pde.gencode = 1;
% generate source codes and store them in app folder
if pde.gencode==1
  %gencode(pde);
  kkgencode(pde);
end

compilerstr = cmakecompile(pde); % use cmake to compile C++ source codes 

runcode(pde, 1); % run C++ code
% steady_PTC(pde, mesh, master, UDG, [], mesh.vdg(:,1,:));
%%
% get solution from output files in dataout folder
sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],2,1);
pde.read_uh = 0;

%%
pde.source = 'source2d';
pde.flux = 'flux_mpp';
pde.fbou = 'fbou_mpp';
pde.fhat = 'fhat_mpp';
pde.arg = {rho_ref, u_ref, rhoE_ref, T_ref, mu_ref, kappa_ref, cp, L_ref, Ec, pde.tau};
pde.bcm  = [2,5,6];  % 2: Wall, 1: Far-field
pde.bcs  = [ui; ui; ui];
pde.bcd  = [1,1,1];  
pde.bcv  = [0; 0; 0];

pde.fc_q = 1;
pde.fc_u = 1;
pde.time = 0;
pde.tdep = 0;
pde.nd = 2;
mesh1 = hdgmesh(mesh, porder);
mesh1.dgnodes(:,3,:) = mesh.vdg(:,1,:);
master = Master(pde);
UDG0 = initu(mesh1,{ui(1),ui(2),ui(3),ui(4),ui(5),ui(6),ui(7),ui(8),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0});
UDG0(:,ns+1,:) = UDG0(:,ns+1,:).*tanh(dist*4);
UDG0(:,ns+2,:) = UDG0(:,ns+2,:).*tanh(dist*4);
% TnearWall = pinf/(gam-1); % Tinf * (Twall/Tref-1) * exp(-10*dist) + Tinf;
rhoe_dim = rhoE_inf - 1/2 * sum(rho_inf) * u_inf^2;
rhoe = rhoe_dim / (rhoE_ref);
UDG0(:,ns+3,:) = rhoe + 1/2 * (UDG0(:,ns+1,:).^2 + UDG0(:,ns+2,:).^2) ./ sum( UDG0(:,1:ns,:), 2 ); 
UH0 = getuhat(UDG0, mesh1.f2t, master.perm, 8);
[UDG,UH] = hdgsolve(master,mesh1,pde,UDG0,UH0,[]);

% %%
% disp("Iter 2")
% mesh.vdg(:,1,:) = 0.04.*tanh(dist*30);
% mesh.udg = sol;
% [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% runcode(pde, 1); % run C++ code
% sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],2,1);
% %%
% %
% disp("Iter 3")
% mesh.vdg(:,1,:) = 0.025.*tanh(dist*30);
% mesh.udg = sol;
% [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% runcode(pde, 1); % run C++ code
% sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],2,1);
% UDG0_3 = real(sol);
% 
% % %%
% % disp("Iter 4")
% % mesh.vdg(:,1,:) = 0.025.*tanh(dist*25);
% % mesh.udg = UDG0_3;
% % [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% % runcode(pde, 1); % run C++ code
% % sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% % figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],2,1);
% % % UDG0_3 = real(sol);
% 
% %%
% disp("Iter 4")
% UDG0_4 = UDG0_3;
% for dd = [5]
%     disp('----------')
%     disp(dd)
%     disp('----------')
%     mesh.vdg(:,1,:) = 0.025.*tanh(dist*dd);
%     mesh.udg = UDG0_4;
%     [pde,mesh,master,dmd] = preprocessing(pde,mesh);
%     runcode(pde, 1); % run C++ code
%     sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
%     figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[]);
%     UDG0_4 = sol;
% end
% %%
% % disp("Iter 3")
% % mesh.vdg(:,1,:) = 0.02.*tanh(dist*2);
% % mesh.udg = sol;
% % [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% % runcode(pde, 1); % run C++ code
% % sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% % figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],2,1);
% 
% % %%
% % disp("Iter 5")
% % mesh.vdg(:,1,:) = 0.03.*tanh(dist*);
% % mesh.udg = UDG_4;
% % [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% % runcode(pde, 1); % run C++ code
% % sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% % figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[]);
% % UDG0_5 = sol;
% % %%
% % % %%
% % % disp("Iter 5")
% % % mesh.vdg(:,1,:) = 0.03.*tanh(dist*5);
% % % pde_tdep.dt = 1e-4*ones(1,2000)*10;   % time step sizes
% % % pde_tdep.soltime = 100:100:length(pde_tdep.dt); % steps at which solution are collected
% % % 
% % % mesh.udg = UDG;
% % % [pde_tdep,mesh,master,dmd] = preprocessing(pde_tdep,mesh);
% % % runcode(pde_tdep, 1); % run C++ code
% % % %%
% % % sol = getsolution('/home/rloekvh/Exasim/build/dataout/out_t1300',dmd,master.npe);
% % % figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[]);
% % % % UDG = sol;
% % 
% % %%
% % mesh.udg = UDG;
% % [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% % runcode(pde, 1); % run C++ code
% % 
% % %%
% % 
% % 
% %%
% mesh1 = hdgmesh(mesh,pde.porder);
% mesh1 = mkcgmesh(mesh1);
% mesh1.dist = dist;
% master_mat = Master(pde);
% for d = 1:mesh.nd+1
% master_mat.shapvl(:,:,d) = master_mat.shapvt(:,:,d)';
% end
% master_mat.gwvl = master_mat.gwe;
% mesh.dist = dist;
% pde.S0=0.2; pde.lambda = 0.04; pde.kappa=3;
% a = avf(mesh1, master_mat, pde, real(sol));
% figure(1); clf; scaplot(mesh, a)
% mesh.vdg(:,1,:) = a;
% mesh.udg = sol;
% [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% runcode(pde, 1); % run C++ code
% 
% % % 
% % % %%
% % % disp("Iter 5")
% % % % mesh.vdg(:,1,:) = /0.025.*tanh(dist*30);.
% % mesh.vdg(:,1,:) = a;
% % mesh.udg = sol;
% % [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% % runcode(pde, 1); % run C++ code
% %%
% [UDG1, UH1] = hdgsolve_avloop(master, master_mat, mesh, mesh1, pde, mesh.dgnodes, sol, [], 0.04, 4);
% %%
% pde.arg = {gam, Minf, 0, Re, Pr, Tref, Twall, pde.tau}; %using matlab fhat, need same pde.arg as matlab
% 
% wid=1;