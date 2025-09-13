% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
pde.platform = "cpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;             % number of MPI processors
pde.hybrid = 1;               % 0 -> LDG, 1 -> HDG
pde.debugmode = 0;
pde.nd = 2;

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;             % polynomial degree
pde.pgauss = 2*pde.porder;
pde.physicsparam = 1;       % unit thermal conductivity
pde.tau = 1.0;              % DG stabilization parameter
pde.linearsolvertol = 1e-8; % GMRES tolerance
pde.ppdegree = 1;          % degree of polynomial preconditioner
pde.RBdim = 0;
pde.GMRESrestart = 50;
pde.saveSolBouFreq = 1;
pde.ibs = 1;

% create a grid of 8 by 8 on the unit square
[mesh.p,mesh.t] = squaremesh(8,8,1,1);
% mesht = mkmesh_square(8,8,pde.porder,0,1,1,1,1);
% mesh.dgnodes = mesht.dgnodes;
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(2,:))<1e-8, @(p) abs(p(1,:)-1)<1e-8, @(p) abs(p(2,:)-1)<1e-8, @(p) abs(p(1,:))<1e-8};
mesh.boundarycondition = [1;1;1;1]; % Set boundary condition for each boundary

% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh,master,dmd,comstr,runstr] = exasim(pde,mesh);

% plot solution
mesh.porder = pde.porder;
mesh.dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);
figure(1); clf; scaplot(mesh,sol(:,1,:),[],2); axis on; axis equal; axis tight;

% [AE, FE, DUDG, DUDG_DUH, D, F, K, H] = uEquationSchurExasim(dataout, pde, npe, npf, nfe);
% [AE0, FE0, DUDG0, DUDG_DUH0, D0, F0, K0, H0] = uEquationSchurExasim(dataout + "0", pde, npe, npf, nfe, 40);
return;

dataout = pde.buildpath + "/dataout/out";
[AE, FE, DUDG, DUDG_DUH] = NewtonExasim(dataout, pde, npe, npf, nfe);
[AE0, FE0, DUDG0, DUDG_DUH0] = NewtonExasim(dataout + "0", pde, npe, npf, nfe);
[AE1, FE1, DUDG1, DUDG_DUH1] = NewtonExasim(dataout + "1", pde, npe, npf, nfe);

i0 = [41:64 33:40 25:32];
i1 = [1:40];
a = squeeze(AE);
a0 = squeeze(AE0);
a1 = squeeze(AE1);
e0=a(:,:,i0(1:32))-a0;
e1=a(:,:,i1(1:32))-a1;
[max(abs(e0(:))) max(abs(e1(:)))]

i0 = [41:64 33:40 25:32];
i1 = [1:40];
a = squeeze(DUDG_DUH);
a0 = squeeze(DUDG_DUH0);
a1 = squeeze(DUDG_DUH1);
e0=a(:,:,i0(1:32))-a0;
e1=a(:,:,i1(1:32))-a1;
[max(abs(e0(:))) max(abs(e1(:)))]

a = squeeze(FE);
a0 = squeeze(FE0);
a1 = squeeze(FE1);
e0=a(:,i0(1:32))-a0;
e1=a(:,i1(1:32))-a1;
[max(abs(e0(:))) max(abs(e1(:)))]

a = squeeze(DUDG);
a0 = squeeze(DUDG0);
a1 = squeeze(DUDG1);
e0=a(:,i0(1:32))-a0;
e1=a(:,i1(1:32))-a1;
[max(abs(e0(:))) max(abs(e1(:)))]

R = assembleRhsExasim(dataout, pde, npf, 144);
R0 = assembleRhsExasim(dataout + "0", pde, npf, 76);
R1 = assembleRhsExasim(dataout + "1", pde, npf, 76);

R = reshape(R, [npf 144]);
R0 = reshape(R0, [npf 76]);
R1 = reshape(R1, [npf 76]);

% visualize the numerical solution of the PDE model using Paraview
% pde.visscalars = {"temperature", 1};  % list of scalar fields for visualization
% pde.visvectors = {"temperature gradient", [2 3]}; % list of vector fields for visualization
% xdg = vis(sol,pde,mesh); % visualize the numerical solution
% disp("Done!");

% master1 = Master(pde);
% mesh1 = hdgmesh(mesh, pde.porder);
% x = mesh1.dgnodes(:,1,:);
% y = mesh1.dgnodes(:,2,:);
% u = x.*x + y.*y;
% UDG = u;
% UDG(:,2,:) = -2*x;
% UDG(:,3,:) = -2*y;
% UDG = 0*UDG;
% UH = inituhat(master1,mesh1.elcon,UDG,1);
% pde.bcm = [1; 1; 1; 1];
% pde.bcs = [1; 1; 1; 1]*0;
% pde.debugmode=1;
% [UDG,UH]= hdgsolve(master1,mesh1,pde,UDG,UH,[]);
% compareexasim(master1, mesh1, pde);


% nge = 25;
% [Mass, Minv, C, E, udg] = qEquationExasim(pde, mesh);
% [Ru, B, D, fg, fg_udg, sg, sg_udg] = uEquationElemExasim(pde, mesh, nge);
% 
% run('/Users/ngoccuongnguyen/Dropbox (MIT)/HDGv1.0/setup.m');
% ngrid = 3;
% mesh1   = mkmesh_square(ngrid+1,ngrid+1,pde.porder,0,1,1,1,1);
% master = mkmaster(mesh1,2*pde.porder);
% app1.hybrid = "hdg";
% app1.bcm = [1; 1; 1; 1];
% app1.bcs = [1; 1; 1; 1]*0;
% app1.flg_q = 1;
% app1.flg_p = 1;
% app1.appname = "poisson";
% app1.dt = [];
% app1.nc = 3;
% app1.ncu = 1;
% app1.nch = 1;
% app1.time = 0;
% app1.tdep = 0;
% app1.wave=0;
% app1.alag=0;
% app1.adjoint=0;
% app1.linearproblem=0;
% app1.fc_u = 1;
% app1.fc_q = 1;
% app.denseblock=0;
% app1.arg = {};
% app1.flux = "flux";
% app1.source = "source";
% app1.fbou   = "fbou";
% app1.fhat   = "fhat";
% [master,mesh1] = preprocess(master,mesh1,app1.hybrid);
% [Mass1, Minv1, C1, E1] = qequation(master, mesh1);
% 
% x = mesh1.dgnodes(:,1,:);
% y = mesh1.dgnodes(:,2,:);
% u = x.*x + y.*y;
% UDG = u;
% UDG(:,2,:) = -2*x;
% UDG(:,3,:) = -2*y;
% SDG = [];
% 
% [Ru1, B1, D1] = uequationelem(master,app1,mesh1.dgnodes,UDG,SDG);
% 
% max(abs(Minv(:)-Minv1(:)))
% max(abs(C(:)-C1(:)))
% max(abs(E(:)-E1(:)))
% max(abs(udg(:)-UDG(:)))
% 
% max(abs(Ru(:)-Ru1(:)))
% max(abs(B(:)-B1(:)))
% max(abs(D(:)-D1(:)))
% 
% UH = inituhat(master,mesh1.elcon,UDG,1);
% [Ru1, Rh1, B1, D1, F1, G1, K1, H1] = uequationelemface(master,app1,mesh1,UDG,UH,SDG);
% [Ru, Rh, B, D, F, G, K, H] = uEquationFaceExasim(pde, mesh);
% 
% max(abs(Ru(:)-Ru1(:)))
% max(abs(Rh(:)-Rh1(:)))
% max(abs(B(:)-B1(:)))
% max(abs(D(:)-D1(:)))
% max(abs(F(:)-F1(:)))
% max(abs(G(:)-G1(:)))
% max(abs(K(:)-K1(:)))
% max(abs(H(:)-H1(:)))
% 
% clear app1;
% app1.hybrid = "hdg";
% app1.bcm = [1; 1; 1; 1];
% app1.bcs = [1; 1; 1; 1]*0;
% app1.flg_q = 1;
% app1.flg_p = 1;
% app1.time = 0;
% app1.tdep = 0;
% app1.wave=0;
% app1.fc_u = 1;
% app1.fc_q = 1;
% app1.denseblock=0;
% app1.arg = {};
% app1.flux = "flux";
% app1.source = "source";
% app1.fbou   = "fbou";
% app1.fhat   = "fhat";
% [master,mesh1] = preprocess(master,mesh1,app1.hybrid);
% 
% pde = initializepde("backend");
% pde.porder = 4;             % polynomial degree
% pde.physicsparam = 1;       % unit thermal conductivity
% pde.pgauss = 2*pde.porder;
% pde.nd = 2;
% pde.arg = {};
% pde.bcm = [1; 1; 1; 1];
% pde.bcs = [1; 1; 1; 1]*0;
% pde.denseblock=0;
% pde.flg_p = 0;
% pde.localsolve = 1;
% pde.adjoint = 0;
% pde.flux = "flux";
% pde.source = "source";
% pde.fbou   = "fbou";
% pde.fhat   = "fhat";
% pde.nc = 3;
% 
% master1 = Master(pde);
% mesh1 = hdgmesh(mesh, pde.porder);
% UH = inituhat(master1,mesh1.elcon,UDG,1);
% [UDG,UH]= hdgsolve(master1,mesh1,pde,UDG,UH,[]);
% 
% 
