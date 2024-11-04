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
pde.platform = "cpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors
pde.hybrid = 1;
% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 4;          % polynomial degree
% pde.torder = 1;          % time-stepping order of accuracy
% pde.nstage = 1;          % time-stepping number of stages
% pde.dt = 1e-4*ones(1,10000)*5;   % time step sizes
% pde.visdt = pde.dt(1);         % visualization timestep size
% pde.saveSolFreq = 100;          % solution is saved every 10 time steps
% pde.soltime = 100:100:length(pde.dt); % steps at which solution are collected

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
pde.read_uh = 0;
% pde.GMRESrestart = 2000;         %try 50
% pde.linearsolvertol = 1e-5; % GMRES tolerance
% pde.linearsolveriter = 50; %try 100
% pde.RBdim = 5;

% pde.ppdegree = 50;
% pde.precMatrixType=2;           % preconditioning type
% pde.ptcMatrixType=0;
% pde.NLtol = 1e-8;              % Newton tolerance
% pde.NLiter = 30;                 % Newton iterations
% pde.matvectol=1e-6;             % tolerance for matrix-vector multiplication
% pde.timestepOffset=0;
% pde.AV = 1;
pde.GMRESrestart = 100;         %try 50
pde.linearsolvertol = 1e-8; % GMRES tolerance
pde.linearsolveriter = 100; %try 100
pde.RBdim = 0;
pde.ppdegree = 40;
pde.precMatrixType=2;           % preconditioning type
pde.ptcMatrixType=0;
pde.NLtol = 1e-6;              % Newton tolerance
pde.NLiter = 30;                 % Newton iterations
pde.matvectol=1e-6;             % tolerance for matrix-vector multiplication
pde.timestepOffset=0;
pde.AV = 1;
%
pde_tdep = pde;
pde_tdep.linearsolveriter = 40;        % number of GMRES iterations
pde_tdep.GMRESrestart = 20;
pde_tdep.linearsolvertol = 1e-3; % GMRES tolerance
pde_tdep.linearsolveriter = 40;
pde_tdep.NLtol = 1e-6;              % Newton tolerance
pde_tdep.NLiter = 3;                 % Newton iterations
pde_tdep.dt = 1e-4*ones(1,1000)*5*10;   % time step sizes
pde_tdep.soltime = 100:100:length(pde_tdep.dt); % steps at which solution are collected
pde_tdep.torder = 1;          % time-stepping order of accuracy
pde_tdep.nstage = 1;          % time-stepping number of stages
pde_tdep.visdt = pde_tdep.dt(1);         % visualization timestep size
pde_tdep.saveSolFreq = 100;          % solution is saved every 10 time steps
%

% [mesh.p,mesh.t,mesh.dgnodes] = mkmesh_circincirc_Ma17b(pde.porder,201,201,1,3,4);
mesh = mkmesh_square(31,21,pde.porder,1,1,1,1,1);
mesh.p(1,:) = logdec(mesh.p(1,:), 6);
mesh.dgnodes(:,1,:) = logdec(mesh.dgnodes(:,1,:), 6);
mesh = mkmesh_halfcircle(mesh, 1, 3, 4.7, pi/2, 3*pi/2);
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

pde.gencode = 0;
% generate source codes and store them in app folder
if pde.gencode==1
  %gencode(pde);
  kkgencode(pde);
end

compilerstr = cmakecompile(pde); % use cmake to compile C++ source codes 

runcode(pde, 1); % run C++ code
%%
% get solution from output files in dataout folder
sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],2,1);
pde.read_uh = 1;

disp("Iter 2")
mesh.vdg(:,1,:) = 0.04.*tanh(dist*30);
mesh.udg = sol;
[pde,mesh,master,dmd] = preprocessing(pde,mesh);
runcode(pde, 1); % run C++ code
sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],2,1);
%%
%
disp("Iter 3")
mesh.vdg(:,1,:) = 0.025.*tanh(dist*30);
mesh.udg = sol;
[pde,mesh,master,dmd] = preprocessing(pde,mesh);
runcode(pde, 1); % run C++ code
sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],2,1);
UDG0_3 = real(sol);

% %%
% disp("Iter 4")
% mesh.vdg(:,1,:) = 0.025.*tanh(dist*25);
% mesh.udg = UDG0_3;
% [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% runcode(pde, 1); % run C++ code
% sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],2,1);
% % UDG0_3 = real(sol);

%%
disp("Iter 4")
UDG0_4 = UDG0_3;
for dd = [5]
    disp('----------')
    disp(dd)
    disp('----------')
    mesh.vdg(:,1,:) = 0.025.*tanh(dist*dd);
    mesh.udg = UDG0_4;
    [pde,mesh,master,dmd] = preprocessing(pde,mesh);
    runcode(pde, 1); % run C++ code
    sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
    figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[]);
    UDG0_4 = sol;
end
%%
% disp("Iter 3")
% mesh.vdg(:,1,:) = 0.02.*tanh(dist*2);
% mesh.udg = sol;
% [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% runcode(pde, 1); % run C++ code
% sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],2,1);

% %%
% disp("Iter 5")
% mesh.vdg(:,1,:) = 0.03.*tanh(dist*);
% mesh.udg = UDG_4;
% [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% runcode(pde, 1); % run C++ code
% sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[]);
% UDG0_5 = sol;
% %%
% % %%
% % disp("Iter 5")
% % mesh.vdg(:,1,:) = 0.03.*tanh(dist*5);
% % pde_tdep.dt = 1e-4*ones(1,2000)*10;   % time step sizes
% % pde_tdep.soltime = 100:100:length(pde_tdep.dt); % steps at which solution are collected
% % 
% % mesh.udg = UDG;
% % [pde_tdep,mesh,master,dmd] = preprocessing(pde_tdep,mesh);
% % runcode(pde_tdep, 1); % run C++ code
% % %%
% % sol = getsolution('/home/rloekvh/Exasim/build/dataout/out_t1300',dmd,master.npe);
% % figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[]);
% % % UDG = sol;
% 
% %%
% mesh.udg = UDG;
% [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% runcode(pde, 1); % run C++ code
% 
% %%
% 
% 
%%
mesh1 = hdgmesh(mesh,pde.porder);
mesh1 = mkcgmesh(mesh1);
mesh1.dist = dist;
master_mat = Master(pde);
for d = 1:mesh.nd+1
master_mat.shapvl(:,:,d) = master_mat.shapvt(:,:,d)';
end
master_mat.gwvl = master_mat.gwe;
mesh.dist = dist;
pde.S0=0.2; pde.lambda = 0.04; pde.kappa=3;
a = avf(mesh1, master_mat, pde, real(sol));
figure(1); clf; scaplot(mesh, a)
mesh.vdg(:,1,:) = a;
mesh.udg = sol;
[pde,mesh,master,dmd] = preprocessing(pde,mesh);
runcode(pde, 1); % run C++ code

% % 
% % %%
% % disp("Iter 5")
% % % mesh.vdg(:,1,:) = /0.025.*tanh(dist*30);.
% mesh.vdg(:,1,:) = a;
% mesh.udg = sol;
% [pde,mesh,master,dmd] = preprocessing(pde,mesh);
% runcode(pde, 1); % run C++ code
%%
[UDG1, UH1] = hdgsolve_avloop(master, master_mat, mesh, mesh1, pde, mesh.dgnodes, sol, [], 0.04, 4);
%%
pde.arg = {gam, Minf, 0, Re, Pr, Tref, Twall, pde.tau}; %using matlab fhat, need same pde.arg as matlab

wid=1;
gam1=gam-1;
Tinf = pinf/(gam-1);
Ttinf = Tinf * (1 + (gam-1)/2 * Minf^2);
TisoW = Twall/Tref * Tinf;
deltaT = Ttinf - TisoW;
elemAvg = 0;
Cp1 = cell(7,1);
Cf1 = cell(7,1);
Ch1 = cell(7,1);
x1 = cell(7,1);
theta1 = cell(7,1);
for i = 1:7
  [Cp1{i},Cf1{i},x1{i},~,~,~,Ch1{i}]=getsurfacedata2_hdgcode(master,mesh1,pde,UDG1{i},UH1{i},wid,elemAvg,deltaT);
  theta1{i} = atan2(x1{i}(:,2),-x1{i}(:,1));
  ii = theta1{i}>0; Cf1{i}(ii) = -Cf1{i}(ii);  
end

ExpCp=[-79.89803449878352, 0.1260397830018084
-51.08143571799568, 0.4360019727108335
-23.173505371641017, 0.8803550879500246
-16.52005139279954, 0.9347197106690777
-7.263606790410322, 1.015403583758014
-4.895984254120992, 0.9964162419858622
-2.651376397583445, 1.0018576360348512
-0.068068122795987, 1.0024823277987833
2.3315382302288072, 1.0012987012987011
4.569585303846259, 1.04670392898241
7.04546075831716, 1.031530494821634
9.424564664716657, 0.9969587374650665
19.536371340313277, 0.884070360019727
26.368606653727344, 0.7807825086306099
33.379623301714005, 0.686059510110143
40.85727563489243, 0.5671214861088278
47.875673163664196, 0.4568633897747822
62.33727891528389, 0.2556962025316456];

ExpCh = [-39.48338121721711, 0.5721278800924818
-36.33780542449722, 0.6621063541417523
-33.402567279053216, 0.6916527146615641
-30.272015603995886, 0.7323447341146456
-27.242679037402144, 0.7949772781631188
-24.13584965339097, 0.8322410906481702
-21.227496771132607, 0.8350155465199712
-14.874667228972823, 0.9798453320577213
0, 1.0031411942916366
-4.832230686101373, 0.9652714661564218
1.3181686391312368, 0.9665151877541257
3.942644771870634, 0.8834250179382921
13.320066422414929, 0.8452204416806186
16.361264134531748, 0.8857689547955034
25.390758849732464, 0.7716176353344495
28.526054982999025, 0.72126285577613
34.67250059305727, 0.6612293709638841
37.803843011149475, 0.6634298014828988
40.72801075410527, 0.5387546838874272
43.93368301747542, 0.5104839352626963
46.887898995756345, 0.5009009009009009
49.92118927752445, 0.42978553775013956
53.01694825904742, 0.3975285019532807
-0.8642821371148423, 0.9632942677190464
-2.3129233769999207, 0.9927290121980387
-5.58659954136904, 0.9211512397353104
-11.97026806188882, 0.8921310691222195
-8.991539049526871, 0.966276010523798
-12.006642241492923, 0.9450689627680777
-9.717441155539154, 0.9076297536474527
-7.160178180763857, 0.9436976799808657];

figure(1); clf; hold on;
plot(theta1{1}*180/pi,-Cp1{1}/max(-Cp1{1}(:)),':','Color',  [0 0.4470 0.7410],'LineWidth',2,'MarkerSize',8);
plot(theta1{3}*180/pi,-Cp1{3}/max(-Cp1{3}(:)),'-.', 'Color', [0.6350 0.0780 0.1840],'LineWidth',2,'MarkerSize',8);
plot(theta1{5}*180/pi,-Cp1{5}/max(-Cp1{5}(:)),'--','Color', [0.4660 0.6740 0.1880], 'LineWidth',2,'MarkerSize',8);
plot(theta1{7}*180/pi,-Cp1{7}/max(-Cp1{7}(:)),'-','Color', [0.4940 0.1840 0.5560],'LineWidth',2,'MarkerSize',6);
plot(ExpCp(:,1),ExpCp(:,2),'ko','LineWidth',2,'MarkerSize',8);
set(gca,'FontSize',18); 
axis tight; axis on; box on; 
xlabel("$\theta$", 'interpreter', 'latex', 'FontSize', 20);
ylabel("$p/p_0$", 'interpreter', 'latex', 'FontSize', 20);
leg = legend(["$n=1$","$n=3$","$n=5$", "$n=7$", "$\mbox{Experiment}$"], 'interpreter', 'latex', 'FontSize', 16, 'Location', 'NW');
leg.ItemTokenSize = [30,10];
axis([-90 90 0 1.05]);
ax = gca;
xticks([-90:30:90]);
yticks([0:0.2:1.0]);
exportgraphics(ax,"nscyl8_adaptive1_pressurecoefficient.png",'Resolution',200); 

figure(2); clf; hold on;
plot(theta1{1}*180/pi,Ch1{1}/max(Ch1{1}(:)),':','Color',  [0 0.4470 0.7410],'LineWidth',2,'MarkerSize',8);
plot(theta1{3}*180/pi,Ch1{3}/max(Ch1{3}(:)),'-.', 'Color', [0.6350 0.0780 0.1840],'LineWidth',2,'MarkerSize',8);
plot(theta1{5}*180/pi,Ch1{5}/max(Ch1{5}(:)),'--','Color', [0.4660 0.6740 0.1880], 'LineWidth',2,'MarkerSize',8);
plot(theta1{7}*180/pi,Ch1{7}/max(Ch1{7}(:)),'-','Color', [0.4940 0.1840 0.5560],'LineWidth',2,'MarkerSize',6);
plot(ExpCh(:,1),ExpCh(:,2),'ko','LineWidth',2,'MarkerSize',8);
set(gca,'FontSize',18); 
axis tight; axis on; box on; 
xlabel("$\theta$", 'interpreter', 'latex', 'FontSize', 20);
ylabel("$q/q_0$", 'interpreter', 'latex', 'FontSize', 20);
leg = legend(["$n=1$","$n=3$","$n=5$", "$n=7$", "$\mbox{Experiment}$"], 'interpreter', 'latex', 'FontSize', 16, 'Location', 'NW');
leg.ItemTokenSize = [30,10];
axis([-90 90 0 1.05]);
ax = gca;
xticks([-90:30:90]);
yticks([0:0.2:1.0]);
exportgraphics(ax,"nscyl8_adaptive1_heatflux.png",'Resolution',200); 


% wid=1;
% gam1=gam-1;
% Tinf = pinf/(gam-1);
% Ttinf = Tinf * (1 + (gam-1)/2 * Minf^2);
% TisoW = Twall/Tref * Tinf;
% deltaT = Ttinf - TisoW;
% elemAvg = 0;
% Cp1 = cell(7,1);
% Cf1 = cell(7,1);
% Ch1 = cell(7,1);
% x1 = cell(7,1);
% theta1 = cell(7,1);
% for i = 1:7
%   [Cp1{i},Cf1{i},x1{i},~,~,~,Ch1{i}]=getsurfacedata2(master,mesh1,UDG1{i}(:,1:pde.ncu,:),pde.physicsparam,wid,elemAvg,deltaT,elemAvg);
%   theta1{i} = atand(x1{i}(:,2)./x1{i}(:,1));
%   ii = theta1{i}>0; Cf1{i}(ii) = -Cf1{i}(ii);  
% end

% ExpCp=[-79.89803449878352, 0.1260397830018084
% -51.08143571799568, 0.4360019727108335
% -23.173505371641017, 0.8803550879500246
% -16.52005139279954, 0.9347197106690777
% -7.263606790410322, 1.015403583758014
% -4.895984254120992, 0.9964162419858622
% -2.651376397583445, 1.0018576360348512
% -0.068068122795987, 1.0024823277987833
% 2.3315382302288072, 1.0012987012987011
% 4.569585303846259, 1.04670392898241
% 7.04546075831716, 1.031530494821634
% 9.424564664716657, 0.9969587374650665
% 19.536371340313277, 0.884070360019727
% 26.368606653727344, 0.7807825086306099
% 33.379623301714005, 0.686059510110143
% 40.85727563489243, 0.5671214861088278
% 47.875673163664196, 0.4568633897747822
% 62.33727891528389, 0.2556962025316456];

% ExpCh = [-39.48338121721711, 0.5721278800924818
% -36.33780542449722, 0.6621063541417523
% -33.402567279053216, 0.6916527146615641
% -30.272015603995886, 0.7323447341146456
% -27.242679037402144, 0.7949772781631188
% -24.13584965339097, 0.8322410906481702
% -21.227496771132607, 0.8350155465199712
% -14.874667228972823, 0.9798453320577213
% 0, 1.0031411942916366
% -4.832230686101373, 0.9652714661564218
% 1.3181686391312368, 0.9665151877541257
% 3.942644771870634, 0.8834250179382921
% 13.320066422414929, 0.8452204416806186
% 16.361264134531748, 0.8857689547955034
% 25.390758849732464, 0.7716176353344495
% 28.526054982999025, 0.72126285577613
% 34.67250059305727, 0.6612293709638841
% 37.803843011149475, 0.6634298014828988
% 40.72801075410527, 0.5387546838874272
% 43.93368301747542, 0.5104839352626963
% 46.887898995756345, 0.5009009009009009
% 49.92118927752445, 0.42978553775013956
% 53.01694825904742, 0.3975285019532807
% -0.8642821371148423, 0.9632942677190464
% -2.3129233769999207, 0.9927290121980387
% -5.58659954136904, 0.9211512397353104
% -11.97026806188882, 0.8921310691222195
% -8.991539049526871, 0.966276010523798
% -12.006642241492923, 0.9450689627680777
% -9.717441155539154, 0.9076297536474527
% -7.160178180763857, 0.9436976799808657];

% figure(1); clf; hold on;
% plot(theta1{1},-Cp1{1}/max(-Cp1{1}(:)),':','Color',  [0 0.4470 0.7410],'LineWidth',2,'MarkerSize',8);
% plot(theta1{3},-Cp1{3}/max(-Cp1{3}(:)),'-.', 'Color', [0.6350 0.0780 0.1840],'LineWidth',2,'MarkerSize',8);
% plot(theta1{5},-Cp1{5}/max(-Cp1{5}(:)),'--','Color', [0.4660 0.6740 0.1880], 'LineWidth',2,'MarkerSize',8);
% plot(theta1{7},-Cp1{7}/max(-Cp1{7}(:)),'-','Color', [0.4940 0.1840 0.5560],'LineWidth',2,'MarkerSize',6);
% plot(ExpCp(:,1),ExpCp(:,2),'ko','LineWidth',2,'MarkerSize',8);
% set(gca,'FontSize',18); 
% axis tight; axis on; box on; 
% xlabel("$\theta$", 'interpreter', 'latex', 'FontSize', 20);
% ylabel("$p/p_0$", 'interpreter', 'latex', 'FontSize', 20);
% leg = legend(["$n=1$","$n=3$","$n=5$", "$n=7$", "$\mbox{Experiment}$"], 'interpreter', 'latex', 'FontSize', 16, 'Location', 'NW');
% leg.ItemTokenSize = [30,10];
% axis([-90 90 0 1.05]);
% ax = gca;
% xticks([-90:30:90]);
% yticks([0:0.2:1.0]);
% exportgraphics(ax,"nscyl8_adaptive1_pressurecoefficient.png",'Resolution',200); 

% figure(2); clf; hold on;
% plot(theta1{1}*180/pi,Ch1{1}/max(Ch1{1}(:)),':','Color',  [0 0.4470 0.7410],'LineWidth',2,'MarkerSize',8);
% plot(theta1{3}*180/pi,Ch1{3}/max(Ch1{3}(:)),'-.', 'Color', [0.6350 0.0780 0.1840],'LineWidth',2,'MarkerSize',8);
% plot(theta1{5}*180/pi,Ch1{5}/max(Ch1{5}(:)),'--','Color', [0.4660 0.6740 0.1880], 'LineWidth',2,'MarkerSize',8);
% plot(theta1{7}*180/pi,Ch1{7}/max(Ch1{7}(:)),'-','Color', [0.4940 0.1840 0.5560],'LineWidth',2,'MarkerSize',6);
% plot(ExpCh(:,1),ExpCh(:,2),'ko','LineWidth',2,'MarkerSize',8);
% set(gca,'FontSize',18); 
% axis tight; axis on; box on; 
% xlabel("$\theta$", 'interpreter', 'latex', 'FontSize', 20);
% ylabel("$q/q_0$", 'interpreter', 'latex', 'FontSize', 20);
% leg = legend(["$n=1$","$n=3$","$n=5$", "$n=7$", "$\mbox{Experiment}$"], 'interpreter', 'latex', 'FontSize', 16, 'Location', 'NW');
% leg.ItemTokenSize = [30,10];
% axis([-90 90 0 1.05]);
% ax = gca;
% xticks([-90:30:90]);
% yticks([0:0.2:1.0]);
% exportgraphics(ax,"nscyl8_adaptive1_heatflux.png",'Resolution',200); 


% mesh.udg = sol;
% 
% [pde_tdep,mesh,master,dmd] = preprocessing(pde_tdep,mesh);
% runcode(pde_tdep, 1); % run C++ code
% sol = fetchsolution(pde,master,dmd, pde.buildpath + '/dataout');
% figure(1); clf; scaplot(mesh, eulereval(sol, 'M',gam,Minf),[],2,1);
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
