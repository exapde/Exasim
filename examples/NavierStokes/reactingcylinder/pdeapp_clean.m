% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
ii = ii(end);
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();
pde.buildpath=string(pwd()); 
pde.exasimpath = string(pwd());
addpath(char(srcdir + "/Modeling/reactingflow/kineticsMatlab/"));
addpath(char(srcdir + "/Modeling/reactingflow/transportMatlab/"));
addpath(char(srcdir + "/Modeling/reactingflow/thermodynamicsMatlab/"));
addpath(char(srcdir + "/Modeling/five_species_air/"));

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors
pde.hybrid = 1;
pmin = 2;
pde.gencode = 1;
N = 80;
M = 40;
mesh_alph = 8;
dfactor = 10;
xc=0.95;

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = pmin;          % polynomial degree
porder = pde.porder;
ns=5;

% % These values are found with a call to Mutation using equilibrium composition at the free-stream
rho_inf = [1.0055500000000001e-09 0.00035318161606 1.5872374699999997e-05 0.00116691299088 1.103201281e-05];
u_inf    = 5956;
T_inf    = 901;
L_ref    = 0.045;
T_wall   = 300;

% % HEG3
[constants] = init_rho_i_T(rho_inf(:), T_inf, u_inf, 0.045);

rho_ref = constants.rho_ref;
u_ref = constants.u_ref;
rhoE_ref = constants.rhoE_ref; 
T_ref = constants.T_ref;
mu_ref = constants.mu_ref;
kappa_ref = constants.kappa_ref;
cp_ref = constants.cp; 

rho_inf = constants.rho_inf;
rhou_inf = constants.rhou_inf;
rhov_inf = constants.rhov_inf;
rhoE_inf = constants.rhoE_inf;

% Reference quantities for nondimensionalizing
% rho_ref = 0.001547;  % sum of inflow densities
avfactor = 10;       % when nondim uses a_ref rather than u_inf, I found I needed to start with a larger viscosity param; perhaps due to conditioning
pde.tau = 10.0;  % We nondimensionalize by the free-stream speed of sound, so tau should be greater than the free-stream Mach

% % Nondimensional constants
Re = rho_ref * L_ref * u_ref / mu_ref;
Pr = mu_ref * cp_ref / kappa_ref;
Ec = u_ref^2 / (cp_ref * T_ref);

% Load into Exasim data structures
pde.physicsparam(1:ns) = rho_inf;
pde.physicsparam(ns+1) = rhou_inf;
pde.physicsparam(ns+2) = rhov_inf;
pde.physicsparam(ns+3) = rhoE_inf;
pde.physicsparam(9)  = Pr;
pde.physicsparam(10) = Re;
pde.physicsparam(11) = Ec;

% Nondimensional terms
pde.externalparam = zeros(1,8);
pde.externalparam(1) = rho_ref;  
pde.externalparam(2) = u_ref;    
pde.externalparam(3) = rhoE_ref; 
pde.externalparam(4) = T_ref;
pde.externalparam(5) = mu_ref;
pde.externalparam(6) = kappa_ref;
pde.externalparam(7) = cp_ref; 
pde.externalparam(8) = L_ref;
pde.externalparam(9) = Ec;
pde.externalparam(10) = Pr;
pde.externalparam(11) = Re;
pde.externalparam(12) = T_wall;
pde.externalparam(13:17) = [0.0, 0.0, 0.0, 0.0, 0.0]; % gamma wall

pde.GMRESrestart = 400;         %WAS200
pde.linearsolvertol = 1e-4; % GMRES tolerance
pde.linearsolveriter = 400; %WAS200
pde.RBdim = 0;
pde.ppdegree = 20;
pde.preconditioner=1;
pde.GMRESortho=1;
pde.NLtol = 1e-6;            % Newton tolerance
pde.NLiter = 30;                 % Newton iterations
pde.matvectol=1e-6;             % tolerance for matrix-vector multiplication
pde.timestepOffset=0;
pde.dae_alpha = 0.0;
pde.dae_beta = 0.0;
% pde.AV = 0;%
% pde.AVsmoothingIter = 0;

% Make mesh
mesh = mkmesh_square(N,M,pde.porder,1,1,1,1,1); % I CAN DO p3 100x100 a4.5
mesh.p(1,:) = logdec(mesh.p(1,:),mesh_alph);
mesh.dgnodes(:,1,:) = logdec(mesh.dgnodes(:,1,:),mesh_alph);
% mesh = mkmesh_square2(N,M,porder,1,1,1,1,1,mesh_alph,xc);
% mesh = mkmesh_halfcircle(mesh,1,3,4.7,pi/2,3*pi/2);
mesh = mkmesh_halfcircle(mesh,1,1.4,2.5,pi/2,3*pi/2); %NOTE:MOSTCASESWERE 1.3

mesh.porder = pde.porder;
mesh.boundaryexpr = {@(p) sqrt(p(1,:).^2+p(2,:).^2)<1+1e-6, @(p) p(1,:)>-1e-7, @(p) abs(p(1,:))<20};
mesh.periodicexpr = {};

% iso-thermal wall, supersonic outflow, supersonic inflow
mesh.boundarycondition = [3;2;1]; 

% distance to the wall
pde.pgauss = 2*(pde.porder);
pde.nd = 2;
pde.elemtype = 1;
master = Master(pde);

[~, ~, jac] = volgeom(master.shapent,permute(mesh.dgnodes,[1 3 2]));
hsz = reshape(sqrt(jac),[],1,size(mesh.dgnodes,3));
[~,cgelcon,rowent2elem,colent2elem,~] = mkcgent2dgent(mesh.dgnodes,1e-8);
hh = dg2cg2(max(hsz,0e-5), cgelcon, colent2elem, rowent2elem);
hh = dg2cg2(hh, cgelcon, colent2elem, rowent2elem);
hh = hh ./ pde.porder;

mesh.f = facenumbering(mesh.p,mesh.t,pde.elemtype,mesh.boundaryexpr,mesh.periodicexpr);
f = mkf(mesh.t,mesh.f,2);
dist = meshdist(f,mesh.dgnodes,master.perm,[1]); % distance to the wall

% intial solution
ui = [ rho_inf(:) / rho_ref; rhou_inf / (rho_ref *u_ref); 0.0; rhoE_inf / (rhoE_ref)];
UDG = initu(mesh,{ui(1),ui(2),ui(3),ui(4),ui(5),ui(6),ui(7),ui(8),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0});
UDG(:,ns+1,:) = UDG(:,ns+1,:).*tanh(dist*dfactor);
UDG(:,ns+2,:) = UDG(:,ns+2,:).*tanh(dist*dfactor);
if mesh.boundarycondition(1) == 5
    UDG(:,1:3,:) = UDG(:,1:3,:) .*tanh(dist*dfactor);
end

% force energy to match free-stream and wall temp with velocity going to 0 at the wall
Tsmooth = T_wall + (T_inf - T_wall) .*tanh(dfactor*dist);
[species_thermo_structs, Mw, RU] = thermodynamicsModels();
rho_init_dim = UDG(:,1:5,:)*rho_ref;

for ie = 1:size(Tsmooth,3)
    for ip = 1:size(Tsmooth,1)
        rhoutmp = UDG(ip,6,ie)*rho_ref*u_ref;
        rhovtmp = UDG(ip,7,ie)*rho_ref*u_ref;
        rhotmp = rho_init_dim(ip,:,ie)';
        Ttmp = Tsmooth(ip,1,ie);
        Ptmp = pressure(Ttmp, rhotmp, Mw);
        Xtmp = X_i(rhotmp, Mw);
        etmp = mixtureEnergyMass(Ttmp, Ptmp, Xtmp, Mw, species_thermo_structs);
        rhoetmp = sum(rhotmp) * etmp;
        rhoEtmp = rhoetmp + 0.5 * (rhoutmp^2 + rhovtmp^2) / sum(rhotmp);
        UDG(ip,8,ie) = rhoEtmp / rhoE_ref;
    end
end

% initial guess for wall-temp; unsure if used
WDG = ones(size(mesh.dgnodes,1),1,size(mesh.dgnodes,3));
for i = 1:1
WDG(:,i,:) = T_inf / pde.externalparam(4);
end

% if p == pmin
mesh.udg = UDG;
mesh.wdg = WDG;%
mesh.vdg = zeros(size(mesh.dgnodes,1),2,size(mesh.dgnodes,3));
mesh.vdg(:,1,:) = avfactor*0.008*tanh(dist*dfactor); % For the free-stream solution with this nondim, I typically need a large amount of viscosity
mesh.vdg(:,2,:) = dist;
% generate input files and store them in datain folder
[pde,mesh,master,dmd] = preprocessing(pde,mesh);

%% I use these routines for the python drivers on Lassen...we can discuss if you want to use them. 
runstr = "!cp ./datain/app.bin ./datain/app_steady.bin";
eval(char(runstr));

% generate the unsteady app
init_unsteady_app(pde, mesh);

[pde,mesh,master,dmd] = preprocessing(pde,mesh);

% Open file for writing
fileID = fopen('./datain/sim_info.txt', 'w');

% Write scalars to file (space-separated)
fprintf(fileID, '%d %d %d\n', mesh.ne, master.npe, mesh.nd);

% Close file
fclose(fileID);

% generate source codes and store them in app folder
if pde.gencode==1
    %gencode(pde);
    kkgencode(pde);
    % hdggencode_fbou_only(pde);
    compilerstr = cmakecompile(pde);
end
%%
param = {};
for i = 1:12
    param{i} = pde.externalparam(i);
end

pde.runmode = 10;
pde.dt = 1e-7*ones(500,1);
pde.soltime = 1:length(pde.dt);

avmag = 1.5;

unext = UDG;
tic
failflag  = 0;
for ii = 1:length(avmag)
    disp("==============")
    disp(" Iter : " + string(ii));
    disp("==============")
    mesh.vdg(:,1,:) = avmag(ii)*tanh(dist*dfactor).*max(hh(:));
    figure(1); clf; scaplot(mesh, mesh.vdg(:,1,:));
    uprev = unext;
    U_eval = permute(uprev, [1 3 2]);
        U_eval = reshape(U_eval, [mesh.ne * master.npe pde.nc]);
    
    w = state2mutation_my_temp(U_eval, param);
    T = reshape(w(:,1), [master.npe, 1, mesh.ne]);
    
    % [unext, ~,~, ~, failflag] = steady_PTC(pde, mesh, master, uprev, T, [], mesh.vdg(:,1,:));
    mesh.udg = uprev;
    mesh.wdg = T;
    [pde,mesh,master,dmd] = preprocessing(pde,mesh);
    runcode(pde,1);
    unext = getsolutions(pde.buildpath + "/dataout" + "/outudg",dmd);

    Uplt = unext; 
    a = mesh.vdg(:,1,:);
    % viz_hpc; 
    if any(isnan(unext(:))) || failflag == 1
        disp("NaN detected. Reverting to previous guess...")
        ii = ii-1;
        unext = uprev;
        break;
    end
end
% wnext = getsolution(pde.buildpath + "/dataout" + "/out_wdg",dmd,master.npe);

toc

UDG0 = unext;
WDG0 = WDG;
mesh_mat = hdgmesh(mesh,pde.porder);
mesh_mat = mkcgmesh(mesh_mat);
mesh_mat.ib = [];
mesh_mat.in = 1:size(mesh_mat.p2,2);

mesh_mat.dist = dist;
master_mat = Master(pde);
for d = 1:mesh.nd+1
master_mat.shapvl(:,:,d) = master_mat.shapvt(:,:,d)';
end
master_mat.gwvl = master_mat.gwe;
mesh.dist = dist;

UDG1 = {}; UH1 = {}; WDG1 = {};

pde.S0=0.02; pde.lambda = 1.5; pde.kappa=1.5;
a = avf(mesh_mat, master_mat, pde, unext);
figure(1); clf; scaplot(mesh, a,[],2); 
% %%ƒ
m = 6;
S0 = 0.2;
eta = 0.9; %m = 4;ƒlam
lambda0 = 0.01*avfactor;
kappa0 = 1.5;

lambda = ones(m,1)*lambda0;
for i = 2:m
    lambda(i) = lambda(i-1)*eta;
end
kappa = ones(m,1)*kappa0;
for i = 2:m
    kappa(i) = 1 + (kappa(i-1)-1)*eta;
end

UDG1 = cell(length(lambda),1);
WDG1 = cell(length(lambda),1);
UH1 = cell(length(lambda),1);
ACG1 = cell(length(lambda),1);

%pde.dt = 1e-2*(2.^(0:11));
%pde.runmode = 0;
% pde.max_iter = 1;
for i = 1:length(lambda)

    pde.lambda = lambda(i);
    pde.S0 = S0;
    pde.kappa = kappa(i);
    if i>1
    mesh.vdg(:,1,:) = avf(mesh_mat, master_mat, pde, UDG1{i-1});    
    else
    mesh.vdg(:,1,:) = avf(mesh_mat, master_mat, pde, UDG0);    
    end

    ACG1{i} = mesh.vdg(:,1,:).*hh; 
    
    if i == 1
        mesh.udg = real(UDG0);
        mesh.wdg = real(WDG0);
    else
        mesh.udg = UDG1{i-1};
    end
    figure(1); clf; scaplot(mesh, eulereval(mesh.udg, 'r',1.4,8),[]);
    figure(2); clf; scaplot(mesh, ACG1{i});
    
    
    [pde,mesh,master,dmd] = preprocessing(pde,mesh);
    runcode(pde, 1); % run C++ code
    
    UDG_i = getsolutions(pde.buildpath + "/dataout" + "/outudg",dmd);
    UDG_i = UDG_i(:,:,:,end);
    WDG_i = getsolutions(pde.buildpath + "/dataout" + "/outwdg",dmd);
    WDG_i = WDG_i(:,:,:,end);

    fileID = fopen(pde.buildpath+"/dataout/outuhat_np0.bin",'r');
    UH_i = fread(fileID,'double');
    UH_i = reshape(UH_i(4:end), pde.ncu, mesh.nf*master.npf, []);
    UH_i = UH_i(:,:,end);
    if any(isnan(UDG_i(:))) ||any(isnan(WDG_i(:)))
    disp("NaN detected. Reverting to previous guess...")
    break
    end

    UDG1{i} = real(UDG_i);
    WDG1{i} = real(WDG_i);
    UH1{i} = real(UH_i);
    i_av = i;
end
save("results.mat","UDG1","WDG1","UH1","mesh","mesh_mat","master","pde","i_av","mesh_alph")