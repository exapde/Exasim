% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();
pde.model = "ModelD";  
pde.modelfile = "pdemodel_poi";

% Choose computing platform and set number of processors
pde.platform = "cpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 4;             % number of MPI processors
pde.hybrid = 1;               % 0 -> LDG, 1 -> HDG
pde.debugmode = 0;

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 2;             % polynomial degree
pde.pgauss = 2*pde.porder;
pde.physicsparam = [1 0 1 1];       % unit thermal conductivity
pde.tau = 1.0;              % DG stabilization parameter
pde.linearsolvertol = 1e-8; % GMRES tolerance
pde.ppdegree = 0;          % degree of polynomial preconditioner
pde.linearsolveriter = 400;
pde.GMRESrestart = 200;
pde.preconditioner = 1;
pde.RBdim = 0;

%then extrude it outwards to create the 3d mesh 
extrLength = 0.05;
ne_z = 3;
zz = linspace(0,extrLength,ne_z+1);

mesh2d = mkmesh_flaredplate2d(2);
pde2d.nd = 2;
pde2d.porder = pde.porder;
pde2d.pgauss = pde.pgauss;
pde2d.elemtype = 1;
pde2d.nodetype = 1;
master2d = Master(pde2d);

[~,cgelcon,rowent2elem,colent2elem,~] = mkcgent2dgent(mesh2d.dgnodes,1e-8);
[~, ~, jac] = volgeom(master2d.shapent,permute(mesh2d.dgnodes,[1 3 2]));
jac = reshape(jac,[],1,size(mesh2d.dgnodes,3));
jac = dg2cg2(jac, cgelcon, colent2elem, rowent2elem);
hm2d = sqrt(dg2cg2(jac, cgelcon, colent2elem, rowent2elem));

mesh = extrudemesh(mesh2d,zz);

%assign boundary conditions to 3d mesh 
L = extrLength;
rnose = 0.0002;
load meshparams.mat
mesh.boundaryexpr = { @(p) (p(3,:)< 1e-6),... % Periodic
@(p) (p(3,:)> L - 1e-6),... % Periodic
@(p) (p(2,:)<yO + 1e-6 & p(1,:)<xO + 1e-6),         ... % Slip wall
@(p) (p(1,:)<xL + 1e-6) &  sqrt(p(1,:).^2+ (p(2,:)-.10625).^2)>2*rnose+1e-6,... % free-stream                       ... % 2 Inflow (cap arc)
@(p) (p(1,:)<xL + 1e-6) &  sqrt(p(1,:).^2+ (p(2,:)-.10625).^2)<2*rnose+1e-6,... % nose                         ... % 3 cap
@(p) (p(2,:)<yL + 1e-6 & p(1,:)>xL - 1e-6 & p(1,:)<xC + 1e-6), ... % Flat plate
@(p) (p(2,:)<yE + 1e-6 & p(1,:)>xC - 1e-6 & p(1,:)<xF + 1e-6),         ... % 5 Flared cone (noslip)
@(p) (p(2,:)<yE + 1e-6 & p(1,:)>xF-.1 - 1e-6 & p(1,:)<xE + 1e-6),         ... % 6 Flared cone (slip)
@(p) (p(1,:)>xT - 1e-6),                                       ... % Outflow
@(p) (p(1,:) > -1e10)                                          ... % Freestream
};
mesh.boundarycondition = [3 3 3 2 1 1 1 1 3 2]; %turbulent 
mesh.periodicexpr = {1, @(p) p([1 2],:), 2, @(p) p([1 2],:)};

%[sol,pde,mesh,master,dmd,comstr,runstr] = exasim(pde,mesh);
pde.Cxxpreprocessing = 0;
[pde,mesh,master,dmd] = preprocessing(pde,mesh);
kkgencode(pde);
compilerstr = cmakecompile(pde); % use cmake to compile C++ source codes 
runstr = runcode(pde, 1); % run C++ code
sol = fetchsolution(pde,master,dmd, pde.buildpath + "/dataout");

npe2 = size(mesh2d.dgnodes,1);
ne2 = size(mesh2d.dgnodes,3);
p1 = pde.porder + 1;
nc = size(sol,2);
sol = reshape(sol, [npe2 p1 nc ne2 ne_z]);

sol2d1 = reshape(sol(:,2,:,:,1), [npe2 nc ne2]);
sol2d2 = reshape(sol(:,2,:,:,end), [npe2 nc ne2]);

% plot solution
mesh.porder = pde.porder;
figure(1); clf;
scaplot(mesh2d,sol2d1(:,1,:),[],2);
axis on; axis equal; axis tight; colorbar;

figure(2); clf;
scaplot(mesh2d,sol2d2(:,1,:),[],2);
axis on; axis equal; axis tight;

return;


