function solhm = hmsmoothing(mesh, sensor, kappa, nhms)

if nargin<3 
    kappa = 1e-2;
end
if nargin<4
    nhms = 1;
end

% initialize pde structure and mesh structure
[pde, mesh0] = initializeexasim();
rpath = pde.buildpath + "/radaptivity";
if exist(rpath,'dir') == 0
    mkdir(rpath);
end

% get only necessary fields
mesh0.p = mesh.p;
mesh0.t = mesh.t;
mesh0.boundaryexpr = mesh.boundaryexpr; 
mesh0.boundarycondition = ones(length(mesh.boundarycondition),1); 
mesh0.porder = mesh.porder;
mesh0.elemtype = mesh.elemtype;
mesh0.f = mesh.f;
mesh0.dgnodes = mesh.dgnodes;
mesh = mesh0;

pde.nd = size(mesh.dgnodes, 2);
pde.elemtype = mesh.elemtype;
pde.porder = mesh.porder;      
pde.nodetype = 1;
pde.pgauss = 2*pde.porder;

pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.platform = "cpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;             % number of MPI processors
pde.hybrid = 1;               % 0 -> LDG, 1 -> HDG
pde.debugmode = 0;

% solve Helmholtz equation with homogeneous Neumann boundary condition
pde.modelfile = "pdemodel_helmholtz";    % name of a file defining the PDE model
pde.physicsparam = kappa;  
pde.tau = 1;              % DG stabilization parameter
pde.linearsolvertol = 1e-8; % GMRES tolerance
pde.preconditioner = 1;
pde.GMRESrestart = 200;
%pde.ppdegree = 10;

master = Master(pde);
mesh.xpe = master.xpe;
mesh.telem = master.telem;

[~,cgelcon,rowent2elem,colent2elem] = mkcgent2dgent(mesh.dgnodes,1e-8);
[~, ~, jac] = volgeom(master.shapent,permute(mesh.dgnodes,[1 3 2]));
jac = reshape(jac,[],1,size(mesh.dgnodes,3));
jac = dg2cg2(jac, cgelcon, colent2elem, rowent2elem);
jac = dg2cg2(jac, cgelcon, colent2elem, rowent2elem);

mesh.vdg = reshape(jac, size(mesh.dgnodes(:,1,:)));
mesh.vdg(:,2,:) = sensor;
%mesh.vdg(:,1,:) = 1;

mesh.boundarycondition(:) = 1;

pde.buildpath = rpath + "/helmholtz"; 
[pde,mesh,master,dmd] = preprocessing(pde,mesh);
if exist(pde.buildpath + "/cpuEXASIM", "file") == 0
  kkgencode(pde);
  cmakecompile(pde); 
end        
runcode(pde, 1); % run C++ code
solhm = fetchsolution(pde,master,dmd,pde.buildpath + '/dataout');

% solve Helmholtz equation again if nhms > 1
for i = 2:nhms
    mesh.vdg(:,2,:) = solhm(:,1,:);
    runcode(pde, 1); % run C++ code
    solhm = fetchsolution(pde,master,dmd,pde.buildpath + '/dataout');
end

% % smoothed mesh density and its gradient
% rho = 1 + 10*solhm(:,1,:);
% drhodx = -10*solhm(:,2,:);
% drhody = -10*solhm(:,3,:);
% % rho = dg2cg(rho, cgelcon, cgent2dgent, rowent2elem);
% % drhodx = dg2cg(drhodx, cgelcon, cgent2dgent, rowent2elem);
% % drhody = dg2cg(drhody, cgelcon, cgent2dgent, rowent2elem);

figure(1); clf; scaplot(mesh,sensor(:,1,:),[],1); axis on; axis equal; 
figure(2); clf; scaplot(mesh,mesh.vdg(:,1,:),[],1); axis on; axis equal; 
figure(3); clf; scaplot(mesh,solhm(:,1,:),[],1); axis on; axis equal; 
figure(4); clf; scaplot(mesh,solhm(:,2,:),[],1); axis on; axis equal; 
figure(5); clf; scaplot(mesh,solhm(:,3,:),[],1); axis on; axis equal; 
