function mesh = mkmesh_seawall3d2(porder)

mesh2d = mkmesh_seawall(porder,1);
[wall, land] = lbseawall();

L1 = 4; R = 2; L2 = 8; L3 = 4;
n = [20 20 40 20 20];
z1 = linspace(0, L1, n(1));
z2 = linspace(L1, L1+1, n(2));
z3 = linspace(L1+1, L1+1+L2, n(3));
z4 = linspace(L1+1+L2, L1+1+L2+1, n(4));
z5 = linspace(L1+1+L2+1, L1+1+L2+1+L3, n(5));

%zz = linspace(0,1,21);
% zz = [z1 z2(2:end) z3(2:end) z4(2:end) z5(2:end)];
zz = [z1 z2(2:end) z3(2:end) z4(2:end) z5(2:end)];
% 
% zz = linspace(0, 10, 60);
% mesh = extrudemesh(mesh2d,zz);
% mesh.p = mesh.p([1 3 2],:);
% mesh.dgnodes = mesh.dgnodes(:,[1 3 2],:);
% 
% zmin = min(mesh2d.p(2,:));
% zmax = max(mesh2d.p(2,:));
% ymax = max(zz);
% 
% mesh.boundaryexpr = {@(p) abs(p(2,:)-0)<1e-3, @(p) abs(p(2,:)-ymax)<1e-3, @(p) abs(p(3,:)-zmin)<1e-3, @(p) abs(p(3,:)-zmax)<1e-3, @(p) p(1,:)<1e-3, @(p) abs(p(1,:)-land(1,1))<1e-3, @(p) abs(p(3,:)-land(2,2))<1e-3, @(p) abs(p(3,:))<100};
% mesh.periodicexpr = [];
% mesh.f = facenumbering(mesh.p,mesh.t,1,mesh.boundaryexpr,mesh.periodicexpr);
% 
% plotboundary(mesh);

%zz = linspace(0, 10, 60);
mesh = extrudemesh(mesh2d,zz);

zmin = min(mesh2d.p(2,:));
zmax = max(mesh2d.p(2,:));
ymax = max(zz);

mesh.boundaryexpr = {@(p) abs(p(3,:)-0)<1e-3, @(p) abs(p(3,:)-ymax)<1e-3, @(p) abs(p(2,:)-zmin)<1e-3, @(p) abs(p(2,:)-zmax)<1e-3, @(p) p(1,:)<1e-3, @(p) abs(p(1,:)-land(1,1))<1e-3, @(p) abs(p(2,:)-land(2,2))<1e-3, @(p) abs(p(2,:))<100};
mesh.periodicexpr = [];
[mesh.f, mesh.tprd, mesh.t2t] = facenumbering(mesh.p,mesh.t,1,mesh.boundaryexpr,mesh.periodicexpr);

plotboundary(mesh);
