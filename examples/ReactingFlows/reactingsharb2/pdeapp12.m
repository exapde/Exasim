% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
ii = ii(end);
run(cdir(1:(ii+5)) + "/install/setpath.m");

load solmach12.mat

% initialize pde structure and mesh structure
[pde,~] = initializeexasim();
% pde.buildpath=string(pwd()); 
% pde.exasimpath = string(pwd());
addpath(char(srcdir + "/Modeling/CNS5air/"));

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel_axial";    % name of a file defining the PDE model

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
pde.tau = 10.0;  % We nondimensionalize by the free-stream speed of sound, so tau should be greater than the free-stream Mach (~26 here)
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

pde.dt = 1e-3*(1.25.^(0:31));
pde.nstage = 1;
pde.torder = 1;
pde.saveSolFreq = 1;

L_ref    = 1;
T_wall   = 1400;
rho_phys_inf = 0.004;
v_phys_inf  = 4000;
T_phys_inf  = 250;
p_phys_inf  = 287*rho_phys_inf*T_phys_inf;

[rho_species_inf, rho_phys_inf, rhov_phys_inf, rhoE_phys_inf] = getEquilibriumState(p_phys_inf, T_phys_inf, v_phys_inf);
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

mesh.boundarycondition = [5 1 1 8 2]; % symmetry, inflow, inflow, wall, outflow
%fb = [f_in, f_out, f_iso, f_slip, f_grad, f_noncat, f_cat, f_cat_gam, f_cat_gam_consistent];

mesh.udg = 0*rho_species;
mesh.udg(:,1:5,:) = rho_species/rho_ref;
mesh.udg(:,6:7,:) = rhov_phys/(rho_ref*v_ref);
mesh.udg(:,8,:) = rhoE_phys/rhoe_ref;
mesh.wdg = T_phys/T_ref;
%mesh.vdg = 1.0*mesh.vdg;
dist = meshdist3(mesh.f,mesh.dgnodes,master.perm,4); % distance to the wall
mesh.vdg = 5e-4*tanh(dist*100);
% qdg = gradu(permute(master.shapen(:,:,2:3),[2 1 3]), mesh.dgnodes, -udg);

% load reactingsol.mat
% mesh.udg = udg;
% mesh.wdg = wdg;

% generate input files and store them in datain folder
[pde,mesh,master,dmd] = preprocessing(pde,mesh);

% generate source codes and store them in app folder
kkgencode(pde);
compilerstr = cmakecompile(pde);
runcode(pde,1);

soludg = getsolutions(pde.buildpath + "/dataout/outudg", dmd);
solwdg = getsolutions(pde.buildpath + "/dataout/outwdg", dmd);
udg = soludg(:,:,:,32);
wdg = solwdg(:,:,:,32);

rho = sum(udg(:,1:5,:),2);
for i = 1:5  
  figure(i); clf; scaplot(mesh, udg(:,i,:)./rho(:,1,:), [], 1); colorbar; colormap('jet')
  hold on; plot(0.0287885, 0.0506758, 'o');
end
figure(6); clf; scaplot(mesh, rho(:,1,:), [], 1); colorbar; colormap('jet');
figure(7); clf; scaplot(mesh, udg(:,6,:)./rho(:,1,:), [], 1); colorbar; colormap('jet')
figure(8); clf; scaplot(mesh, udg(:,7,:)./rho(:,1,:), [], 1); colorbar; colormap('jet')
figure(9); clf; scaplot(mesh, udg(:,8,:), [], 1); colorbar; colormap('jet');
figure(10); clf; scaplot(mesh, wdg(:,1,:), [], 1); colorbar; colormap('jet');

r1 = 1./rho;
rx = sum(udg(:,9:13,:),2);
rux = udg(:,14,:);
ry = sum(udg(:,17:21,:),2);
rvy = udg(:,23,:);
u = udg(:,6,:)./rho;
v = udg(:,7,:)./rho;

ux  = (rux - rx.*u).*r1;    
vy  = (rvy - ry.*v).*r1;            
div = ux + vy - v./mesh.dgnodes(:,2,:);    
figure(11); clf; scaplot(mesh, div, [-100 100], 1); colorbar; colormap('jet');

[~,cgelcon,rowent2elem,colent2elem,~] = mkcgent2dgent(mesh.dgnodes,1e-8);
[~, ~, jac] = volgeom(master.shapent,permute(mesh.dgnodes,[1 3 2]));
jac = reshape(jac,[],1,size(mesh.dgnodes,3));
jac = dg2cg2(jac, cgelcon, colent2elem, rowent2elem);
jac = dg2cg2(jac, cgelcon, colent2elem, rowent2elem);
hm = jac.^0.5;

% limit compression
sigm = 100.0;
comp = limiting(div,0,sigm,1e3,0);
DucrosRatio = 1.0;
c_star = 1.0;
sb = (hm./porder).^0.25 .* (comp/c_star) * DucrosRatio;
avField = 6e-4*limiting(sb,0,4,1e3,0.1);

figure(12); clf; scaplot(mesh, avField, [], 1); colorbar; colormap('jet');


