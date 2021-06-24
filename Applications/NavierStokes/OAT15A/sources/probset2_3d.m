pdemodel = 2;
appname = 'ns';
gpu = 0;
mpiprocs = 16;
porder = 2;
nstage = 3;
torder = 3;
dt = 1e-4;
ntime = 100000;
extrLength = 0.1;
ne_z = 2;
gam = 1.4;
Re = 3e6;
Pr = 0.72;
Minf = 0.73;
pinf = 1/(gam*Minf^2);
tau = 4;
alpha = 3.5*pi/180; % Angle between freestream and x-axis
ui = [1.0, cos(alpha), sin(alpha), 0, 0.5+pinf/(gam-1)];
rmin = 0.4; % minimum density
pmin = 0.45; % minimum presssure
avr = 1e-2/3; % nominal viscosity
avs = 1e-4; % nominal viscosity
%periodicexpr = {};
periodicexpr = {1,'p(:,[1 2])',2,'p(:,[1 2])'};
% generate mesh
load('sol2d.mat');
if (porder==3)
    mesh = meshp3;
    UDG = UDGp3;
    avfix = avfixp3;
end
ne2d = mesh.ne;
zz = linspace(0,extrLength,ne_z+1);
mesh = mkmesh_3dextrudemesh2(mesh,zz,mesh.elemtype);
% initial solution
UDG = extrudesol(UDG, mesh.porder, ne_z);
size(UDG)
UDG(:,16:20,:) = 0;
ind1 = [15 13 12 11 10 8 7 6 5 3 2 1];
ind2 = 12:-1:1;
for i = 1:length(ind1)
    UDG(:,ind1(i),:) = UDG(:,ind2(i),:);
end
UDG(:,[4 9 14],:) = 0;
avfix = extrudesol(avfix, mesh.porder, ne_z);
% make app structure
clear app;
app.version = 'exasim4.0';
%app.compilerflags = '-D TIMING';
app.pdemodel = pdemodel;
app.appname = appname;
app.porder = porder;
app.torder = torder;
app.nstage = nstage;
app.mpiprocs = mpiprocs;
app.elemtype = mesh.elemtype;
app.nodetype = mesh.nodetype;
app.physicsparam = [gam Re Pr Minf rmin avr pmin avs];
app.boundarycondition = [2;2;2;1]; % 1: wall, 2: far field
app.uinf = ui;
app.tau = tau;
app.dt = [loginc(linspace(1e-5,dt,21),2)'; ones(50000,1)*dt];
app.dt=ones(ntime,1)*dt;
app.tdfunc=0;
app.source=0;
app.saveSolFreq = 50;
app.saveSolOpt = 0;
app.matvecorder=1;
app.RBdim = 6;
app.PTCiter = 3;
app.PTCtol = 2e-9;
app.PTCparam = 0.0;
app.GMRESrestart=7;
app.linearsolvertol=0.01;
app.linearsolveriter=8;
app.precMatrixType=2;
app.ptcMatrixType=0;
app.AV = 3;
%app.runmode=1;app.dt=0;
npl = size(avfix,1);
ODG = zeros(npl,13,mesh.ne);
ODG(:,4,:) = avfix;
[~,Minv] = meshmetric(mesh.dgnodes(:,:,1:ne2d),mesh.porder,mesh.elemtype,mesh.nodetype);
ODG = reshape(ODG,[npl,13,ne2d,ne_z]);
for i = 1:ne_z
    ODG(:,5:13,:,i) = Minv;
end
ODG = reshape(ODG,[npl,13,ne2d*ne_z]);
% generate binary files
app = preprocessing(app,mesh.p,mesh.t,mesh.dgnodes,UDG,ODG,[],mesh.bndexpr,periodicexpr);
mkdir('data');
movefile('*.bin','data');
