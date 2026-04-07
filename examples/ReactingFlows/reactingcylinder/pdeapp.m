% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
ii = ii(end);
run(cdir(1:(ii+5)) + "/install/setpath.m");

load solp2.mat

% initialize pde structure and mesh structure
[pde,~] = initializeexasim();
% pde.buildpath=string(pwd()); 
% pde.exasimpath = string(pwd());
addpath(char(srcdir + "/Modeling/CNS5air/"));

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
pde.platform = "cpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 4;             % number of MPI processors
pde.hybrid = 1;               % 0 -> LDG, 1 -> HDG
pde.debugmode = 0;
pde.nd = 2;
pde.porder = 2;               % polynomial degree
pde.pgauss = 2*(pde.porder);
pde.nd = 2;
pde.elemtype = 1;
pde.tau = 5.0;  % We nondimensionalize by the free-stream speed of sound, so tau should be greater than the free-stream Mach (~26 here)
pde.gencode = 1;

pde.GMRESrestart = 250;         %try 50
pde.GMRESortho = 1;
pde.linearsolvertol = 1e-6; % GMRES tolerance
pde.linearsolveriter = 500; %try 100
pde.preconditioner = 1;
pde.RBdim = 0;
pde.ppdegree = 0;
pde.NLtol = 1e-6;              % Newton tolerance
pde.NLiter = 10;                % Newton iterations
pde.matvectol=1e-6;             % tolerance for matrix-vector multiplication
pde.dae_alpha = 0;
pde.dae_beta = 0;
pde.dae_gamma = 0;

pde.dt = []; %1e-3*(1.25.^(0:37));
pde.nstage = 1;
pde.torder = 1;
pde.saveSolFreq = 1;

L_ref    = 0.045;
T_wall   = 300;
rho_species_inf = [1.0055500000000001e-09 0.00035318161606 1.5872374699999997e-05 0.00116691299088 1.103201281e-05];
rho_phys_inf = sum(rho_phys_species_inf);
v_phys_inf  = 5956;
T_phys_inf  = 901;
%p_phys_inf  = 287*rho_phys_inf*T_phys_inf;
%[~, ~, rhov_phys_inf, rhoE_phys_inf] = getEquilibriumState(p_phys_inf, T_phys_inf, v_phys_inf);

rhov_phys_inf = v_phys_inf*rho_phys_inf;
[rhoE_phys_inf, p_phys_inf, ~] = energyFromSpecies(rho_species_inf, T_phys_inf, v_phys_inf, 1e4);
[rho_ref, v_ref, rhoe_ref, p_ref, T_ref, mu_ref, kappa_ref, lambda_ref, cp_ref, cv_ref] = getReferenceState(p_phys_inf, T_phys_inf, v_phys_inf);

% Nondimensional constants
Re = rho_ref * L_ref * v_ref / mu_ref;
Pr = mu_ref * cp_ref / kappa_ref;
Ec = v_ref^2 / (cp_ref * T_ref);

U_ref = [rho_ref v_ref rhoe_ref T_ref mu_ref kappa_ref cp_ref L_ref Ec Pr Re T_wall];
U_inf = [rho_species_inf/rho_ref rhov_phys_inf/(rho_ref*v_ref) 0 rhoE_phys_inf/rhoe_ref];
pde.physicsparam = U_ref;
pde.externalparam = U_inf;
pde.externalparam(9:13) = U_inf(1:5);
pde.externalparam(14:18) = 0;

%fb = [f_in, f_out, f_iso, f_slip, f_grad, f_noncat, f_cat, f_cat_gam, f_cat_gam_consistent];
% iso-thermal wall, supersonic outflow, supersonic inflow
mesh = mkmesh_cylns(pde.porder);
mesh.boundarycondition = [8;2;1]; 

u_ref = v_ref;
dfactor = 10;
ns = 5;

% ui = U_inf;
% UDG = initu(mesh,{ui(1),ui(2),ui(3),ui(4),ui(5),ui(6),ui(7),ui(8)});
% UDG(:,ns+1,:) = UDG(:,ns+1,:).*tanh(dist*dfactor);
% UDG(:,ns+2,:) = UDG(:,ns+2,:).*tanh(dist*dfactor);
% 
% % force energy to match free-stream and wall temp with velocity going to 0 at the wall
% Tsmooth = T_wall + (T_phys_inf - T_wall) .*tanh(dist*dfactor);
% [species_thermo_structs, Mw, RU] = thermodynamicsModels();
% rho_init_dim = UDG(:,1:5,:)*rho_ref;
% 
% for ie = 1:size(Tsmooth,3)
%     for ip = 1:size(Tsmooth,1)
%         rhoutmp = UDG(ip,6,ie)*rho_ref*u_ref;
%         rhovtmp = UDG(ip,7,ie)*rho_ref*u_ref;
%         rhotmp = rho_init_dim(ip,:,ie)';
%         Ttmp = Tsmooth(ip,1,ie);
%         Ptmp = pressure(Ttmp, rhotmp, Mw);
%         Xtmp = X_i(rhotmp, Mw);
%         etmp = mixtureEnergyMass(Ttmp, Ptmp, Xtmp, Mw, species_thermo_structs);
%         rhoetmp = sum(rhotmp) * etmp;
%         rhoEtmp = rhoetmp + 0.5 * (rhoutmp^2 + rhovtmp^2) / sum(rhotmp);
%         UDG(ip,8,ie) = rhoEtmp / rhoe_ref;
%     end
% end
% 
% % initial guess for wall-temp; unsure if used
% WDG = ones(size(mesh.dgnodes,1),1,size(mesh.dgnodes,3));
% for i = 1:1
%  WDG(:,i,:) = Tsmooth/T_ref;
% end

mesh.udg = UDG;
mesh.wdg = WDG;%
mesh.vdg = zeros(size(mesh.dgnodes,1),1,size(mesh.dgnodes,3));
mesh.vdg(:,1,:) = 4e-3*tanh(dist*dfactor); % AV for species
mesh.vdg(:,2,:) = 4e-3*tanh(dist*dfactor); 


% generate input files and store them in datain folder
[pde,mesh,master,dmd] = preprocessing(pde,mesh);

% generate source codes and store them in app folder
kkgencode(pde);
compilerstr = cmakecompile(pde);
runcode(pde,1);

UDG = getsolutions(pde.buildpath + "/dataout/outudg", dmd);
WDG = getsolutions(pde.buildpath + "/dataout/outwdg", dmd);

UDG = UDG(:,:,:,end);
WDG = WDG(:,:,:,end);
rho = 0;
for i = 1:5  
  figure(i); clf; scaplot(mesh, UDG(:,i,:), [], 1); colorbar; colormap('jet')
  rho = rho + UDG(:,i,:);
end
figure(6); clf; scaplot(mesh, UDG(:,6,:), [], 1); colorbar; colormap('jet')
figure(7); clf; scaplot(mesh, UDG(:,7,:), [], 1); colorbar; colormap('jet')
figure(8); clf; scaplot(mesh, UDG(:,8,:), [], 1); colorbar; colormap('jet');
figure(9); clf; scaplot(mesh, WDG(:,1,:)*T_ref, [], 1); colorbar; colormap('jet');
figure(10); clf; scaplot(mesh, rho*rho_ref, [], 1); colorbar; colormap('jet');

