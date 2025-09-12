
m      = 21;
n      = 31;
porder = 4;
nstage = 1;
torder = 1;
hybrid = 'hdg';

gam = 1.4;
epslm = 0.1;
Minf = 0.2;                  % Infinity conditions
pinf = 1/(gam*Minf^2);
alpha = 0*pi/180;
Re = 100;
Pr = 0.72;
tau = 3.0;

ntime = 4;
%dt = linspace(0.01,1,ntime);
dt = [1e-4 1e-3 1e-2 1e-1]*1e-1;

ui = [ 1, cos(alpha), sin(alpha), 0.5+pinf/(gam-1)];
%uip = [pinf, pinf, pinf, pinf];

clear app;
app.source = 'source';
app.flux = 'flux';
app.fbou = 'fbou';
app.fhat = 'fhat';

app.getdqdg = 1;
app.denseblock = 0;
app.hybrid = hybrid;
app.localsolve = 1;
app.arg = {gam,epslm,Re,Pr,Minf,tau};
app.bcm  = [2,1];  % 2: Wall, 1: Far-field
app.bcs  = [ui;ui];

app.bcd  = [1,3];  % 2: Wall, 1: Far-field
app.bcv  = [ui;ui];

app.wave = false;
app.tdep = true;
app.alag = false;
app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;
app.fc_q = 1;
app.fc_u = 1;
app.fc_p = 0;

app.np   = 2;
app.nd   = 2;
app.ncu  = 2+app.nd;               % Number of components of U
app.nch  = app.ncu;                % Number of componets of UH
app.nq   = app.ncu*app.nd;         % Number of componets of Q
app.nc   = app.ncu+app.nq;         % Number of componeents of UDG

app.time = [];
app.dtfc = [];
app.alpha = [];

%mesh = mkmesh_trefftz(m,n,porder,[0.05,0.05,1.98]);
mesh = mkmesh_naca0012(porder,1,1);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

UDG0 = initu(mesh,{ui(1),ui(2),ui(3),ui(4); 0,0,0,0; 0,0,0,0});

UDG0(:,2,:) = UDG0(:,2,:).*tanh(meshdist(mesh,1)*5);
UDG0(:,3,:) = UDG0(:,3,:).*tanh(meshdist(mesh,1)*5);
TnearWall = pinf/(gam-1); 
UDG0(:,4,:) = TnearWall + 0.5*(UDG0(:,2,:).*UDG0(:,2,:) + UDG0(:,3,:).*UDG0(:,3,:));

UH0 = inituhat(master,mesh.elcon,UDG0,app.ncu);

fprintf('\nSteady State Solve\n');

app.fc_q = 1;
app.fc_u = 0;
app.tdep = false;
app.adjoint = 0;
[UDG,UH] = hdg_solve(master,mesh,app,UDG0,UH0,[]);

%figure(1); clf; scaplot(mesh,eulereval(UDG(:,1:app.nch,:),'M',gam),[],2); axis off; %axis equal; axis tight;

%mesh.boundaryexpr = {@(p) abs(p(2,:))<1e-8, @(p) abs(p(1,:)-1)<1e-8, @(p) abs(p(2,:)-1)<1e-8, @(p) abs(p(1,:))<1e-8};

% expressions for domain boundaries
% mesh1 = mesh;
mesh1.t = mesh.t';
mesh1.p = mesh.p';
mesh1.boundaryexpr = {@(p) sqrt((p(1,:)-.5).^2+p(2,:).^2)<3, @(p) abs(p(1,:))<20};
mesh1.boundarycondition = [1;2];
% expressions for curved boundaries
mesh1.curvedboundary = [1 0];
mesh1.curvedboundaryexpr = {@(p) p(2,:).^2-(5*0.01*12*(0.29690*sqrt(abs(p(1,:)))-0.12600*p(1,:)-0.35160*p(1,:).^2+0.28430*p(1,:).^3-0.10150*p(1,:).^4)).^2, @(p) 0};
mesh1.periodicexpr =  {};
mesh1 = hdgmesh(mesh1, master.porder);

app.porder = porder;
app.pgauss = 2*porder;
app.elemtype = 1;
app.nodetype = 1;
master1 = Master(app);
mesh1.telem = master1.telem;
mesh1.tface = master1.telem;
mesh1.xpe = master1.xpe;
mesh1.xpf = master1.xpf;

f = mkf(mesh1.t,mesh1.f,master1.nd);
UDG0 = initu(mesh1,{ui(1),ui(2),ui(3),ui(4); 0,0,0,0; 0,0,0,0});
UDG0(:,2,:) = UDG0(:,2,:).*tanh(meshdist(f,mesh1.dgnodes,master.perm,1)*5);
UDG0(:,3,:) = UDG0(:,3,:).*tanh(meshdist(f,mesh1.dgnodes,master.perm,1)*5);
TnearWall = pinf/(gam-1); 
UDG0(:,4,:) = TnearWall + 0.5*(UDG0(:,2,:).*UDG0(:,2,:) + UDG0(:,3,:).*UDG0(:,3,:));
UH0 = inituhat(master,mesh1.elcon,UDG0,app.ncu);

app.denseblock=0;
[UDG,UH] = hdgsolve(master1,mesh1,app,UDG0,UH0,[]);


