function mesh = mkmesh_seawall3d3(porder)

mesh2d = mkmesh_seawall(porder,1);
[wall, land] = lbseawall();

% L1 = 11; R = 4; L2 = 21; L3 = 11;
% n = [40 40 80 40 40]/2;

L1 = 41; R = 4; L2 = 210; L3 = 41;
n = [40 10 200 10 40]*2;

z1 = linspace(0, L1, n(1));
z2 = linspace(L1, L1+1, n(2));
z3 = linspace(L1+1, L1+1+L2, n(3));
z4 = linspace(L1+1+L2, L1+1+L2+1, n(4));
z5 = linspace(L1+1+L2+1, L1+1+L2+1+L3, n(5));

zz = [z1 z2(2:end) z3(2:end) z4(2:end) z5(2:end)] - (L1 + 1);

mesh = extrudemesh(mesh2d,zz);

ymin = min(mesh2d.p(2,:));
ymax = max(mesh2d.p(2,:));
zmin = min(zz);
zmax = max(zz);

mesh.boundaryexpr = {@(p) abs(p(3,:)-zmin)<1e-3, @(p) abs(p(3,:)-zmax)<1e-3, @(p) abs(p(2,:)-ymin)<1e-3, @(p) abs(p(2,:)-ymax)<1e-3, @(p) p(1,:)<1e-3, @(p) abs(p(1,:)-land(1,1))<1e-3, @(p) abs(p(2,:)-land(2,2))<1e-3, @(p) abs(p(2,:))<100};
mesh.periodicexpr = [];
[mesh.f, mesh.tprd, mesh.t2t] = facenumbering(mesh.p,mesh.t,1,mesh.boundaryexpr,mesh.periodicexpr);

plotboundary(mesh);

x = mesh.p(1,:);
z = mesh.p(3,:);
ind1 = ((-L1-1 <= z) & (z <= -1)); 
ind2 = ((-1 <= z) & (z <= 0)); 
ind3 = ((0 <= z) & (z <= L2)); 
ind4 = ((L2 <= z) & (z <= L2+1)); 
ind5 = (L2+1 <= z) & (z <= L2+1+L3); 

[min(z(ind1)) max(z(ind1))]
[min(z(ind2)) max(z(ind2))]
[min(z(ind3)) max(z(ind3))]
[min(z(ind4)) max(z(ind4))]
[min(z(ind5)) max(z(ind5))]

X = mesh.dgnodes(:,1,:);
Z = mesh.dgnodes(:,3,:);
ine1 = (-L1-1 <= Z) & (Z <= -1); 
ine2 = (-1 <= Z) & (Z <= 0); 
ine3 = (0 <= Z) & (Z <= L2); 
ine4 = (L2 <= Z) & (Z <= L2+1); 
ine5 = (L2+1 <= Z) & (Z <= L2+1+L3); 

[min(Z(ine1)) max(Z(ine1))]
[min(Z(ine2)) max(Z(ine2))]
[min(Z(ine3)) max(Z(ine3))]
[min(Z(ine4)) max(Z(ine4))]
[min(Z(ine5)) max(Z(ine5))]

z1 = z(ind1);
x1 = x(ind1);
z2 = z(ind2);
x2 = x(ind2);
z3 = z(ind3);
x3 = x(ind3);
z4 = z(ind4);
x4 = x(ind4);
z5 = z(ind5);
x5 = x(ind5);

Z1 = Z(ine1);
X1 = X(ine1);
Z2 = Z(ine2);
X2 = X(ine2);
Z3 = Z(ine3);
X3 = X(ine3);
Z4 = Z(ine4);
X4 = X(ine4);
Z5 = Z(ine5);
X5 = X(ine5);

zmax = max(z);
xmax = max(x);
% xc = L1;
% yc = ymax + R;

% map z in [L2+1 L2+1+L3] -> x in [xmax+R xmax+R+L3] 
mesh.p(1,ind5) = xmax+R + (z5-(L2+1));
% map x in [0 xmax] -> z in [L2+R+xmax L2+R] 
mesh.p(3,ind5) = L2+R+xmax - x5; 

XT = X; ZT = Z; 
XT(ine5) = xmax+R + (Z5-(L2+1));
ZT(ine5) = L2+R+xmax - X5; 

% map [L2 L2+1] to [pi pi/2] 
t = pi - (pi/2)*(z4 - L2);
rt = xmax+R - x4;
xt = xmax+R + rt.*cos(t);
zt = L2 + rt.*sin(t);
mesh.p(1,ind4) = xt;
mesh.p(3,ind4) = zt;

t = pi - (pi/2)*(Z4 - L2);
rt = xmax+R - X4;
XT(ine4) = xmax+R + rt.*cos(t);
ZT(ine4) = L2 + rt.*sin(t);
% mesh.dgnodes(:,1,:) = XT;
% mesh.dgnodes(:,3,:) = ZT;

% map [-1 0] to [3/2*pi pi] 
t = 3/2*pi - (pi/2)*(z2 + 1);
rt = xmax+R - x2;
xt = xmax+R + rt.*cos(t);
zt = 0 + rt.*sin(t);
mesh.p(1,ind2) = xt;
mesh.p(3,ind2) = zt;

t = 3/2*pi - (pi/2)*(Z2 + 1);
rt = xmax+R - X2;
XT(ine2) = xmax+R + rt.*cos(t);
ZT(ine2) = 0 + rt.*sin(t);

% map z in [-1 -L-1] -> x in [xmax+R xmax+R+L1] 
mesh.p(1,ind1) = xmax+R - (z1+1);
% map x in [0 xmax] -> z in [-R-xmax -R] 
mesh.p(3,ind1) = -R-xmax + x1; 
XT(ine1) = xmax+R - (Z1+1);
ZT(ine1) = -R-xmax + X1; 

mesh.dgnodes(:,1,:) = XT;
mesh.dgnodes(:,3,:) = ZT;

%mesh.p = mesh.p/0.3048;

plotboundary(mesh,3);
% X = mesh.dgnodes(:,1,:);
% Y = mesh.dgnodes(:,2,:);
% Z = mesh.dgnodes(:,3,:);
% hold on;
% plot3(X(:),Y(:),Z(:),'o');