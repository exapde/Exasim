function mesh = mkmesh_seawall3d(porder)

mesh2d = mkmesh_seawall(porder,1);
[wall, land] = lbseawall();

L1 = 4; R = 2; L2 = 8; L3 = 4;
n = [20 20 40 20 20]/2;
z1 = linspace(0, L1, n(1));
z2 = linspace(L1, L1+1, n(2));
z3 = linspace(L1+1, L1+1+L2, n(3));
z4 = linspace(L1+1+L2, L1+1+L2+1, n(4));
z5 = linspace(L1+1+L2+1, L1+1+L2+1+L3, n(5));

%zz = linspace(0,1,21);
% zz = [z1 z2(2:end) z3(2:end) z4(2:end) z5(2:end)];
zz = [z1 z2(2:end) z3(2:end) z4(2:end) z5(2:end)];

mesh = extrudemesh(mesh2d,zz);
mesh.p = mesh.p([1 3 2],:);
mesh.dgnodes = mesh.dgnodes(:,[1 3 2],:);

zmin = min(mesh2d.p(2,:));
zmax = max(mesh2d.p(2,:));
ymax = max(zz);


mesh.boundaryexpr = {@(p) abs(p(2,:)-0)<1e-3, @(p) abs(p(2,:)-ymax)<1e-3, @(p) abs(p(3,:)-zmin)<1e-3, @(p) abs(p(3,:)-zmax)<1e-3, @(p) p(1,:)<1e-3, @(p) abs(p(1,:)-land(1,1))<1e-3, @(p) abs(p(3,:)-land(2,2))<1e-3, @(p) abs(p(3,:))<100};
mesh.boundarycondition = [1;1;3;2]; 
mesh.periodicexpr = [];
[mesh.f, mesh.tprd, mesh.t2t] = facenumbering(mesh.p,mesh.t,1,mesh.boundaryexpr,mesh.periodicexpr);

plotboundary(mesh);

mesh.p = mesh.p([2 1 3],:);
mesh.dgnodes = mesh.dgnodes(:,[2 1 3],:);
% xmax = max(mesh.p(1,:));
% mesh.p(1,:) = mesh.p(1,:) - xmax;
% mesh.dgnodes(:,1,:) = mesh.dgnodes(:,1,:) - xmax;

plotboundary(mesh,2);

x = mesh.p(1,:);
y = mesh.p(2,:);
ind2 = ((L1 <= x) & (x <= L1+1)); 
ind3 = ((L1+1 <= x) & (x <= L1+1+L2)); 
ind4 = ((L1+1+L2 <= x) & (x <= L1+1+L2+1)); 
ind5 = (L1+1+L2+1 <= x) & (x <= L1+1+L2+1+L3); 

[min(x(ind2)) max(x(ind2))]
[min(x(ind3)) max(x(ind3))]
[min(x(ind4)) max(x(ind4))]
[min(x(ind5)) max(x(ind5))]

X = mesh.dgnodes(:,1,:);
Y = mesh.dgnodes(:,2,:);
ine2 = (L1 <= X) & (X <= L1+1); 
ine3 = (L1+1 <= X) & (X <= L1+1+L2); 
ine4 = (L1+1+L2 <= X) & (X <= L1+1+L2+1); 
ine5 = (L1+1+L2+1 <= X) & (X <= L1+1+L2+1+L3); 

[min(X(ine2)) max(X(ine2))]
[min(X(ine3)) max(X(ine3))]
[min(X(ine4)) max(X(ine4))]
[min(X(ine5)) max(X(ine5))]

x2 = x(ind2);
y2 = y(ind2);
x3 = x(ind3);
y3 = y(ind3);
x4 = x(ind4);
y4 = y(ind4);
x5 = x(ind5);
y5 = y(ind5);

X2 = X(ine2);
Y2 = Y(ine2);
X3 = X(ine3);
Y3 = Y(ine3);
X4 = X(ine4);
Y4 = Y(ine4);
X5 = X(ine5);
Y5 = Y(ine5);

xmax = max(x);
ymax = max(y);
xc = L1;
yc = ymax + R;

% map [L1 L1+1] to [3/2*pi 2*pi] 
t = 3/2*pi + (2*pi - 3/2*pi)*(x2 - L1);
% map [L1 L1+1] x [0 ymax] to circular arc
rt = (yc-y2);
xt = xc + rt.*cos(t);
yt = yc + rt.*sin(t);
mesh.p(1,ind2) = xt;
mesh.p(2,ind2) = yt;

t = 3/2*pi + (2*pi - 3/2*pi)*(X2 - L1);
rt = (yc-Y2);
XT = X; XT(ine2) = xc + rt.*cos(t);
YT = Y; YT(ine2) = yc + rt.*sin(t);
mesh.dgnodes(:,1,:) = XT;
mesh.dgnodes(:,2,:) = YT;

% map y in [0 ymax] -> x in [L1+R+ymax L1+R] 
mesh.p(1,ind3) = L1+R+ymax - y3;
% map x in [L1+1 L1+1+L2] -> y in [ymax+R ymax+R+L2] 
mesh.p(2,ind3) = ymax+R + (x3-L1-1); 

XT = X; XT(ine3) = L1+R+ymax - Y3;
YT = Y; YT(ine3) = ymax+R + (X3-L1-1); 
mesh.dgnodes(:,1,:) = XT;
mesh.dgnodes(:,2,:) = YT;

% map [L1+1+L2 L1+1+L2+1] to [0 pi/2] 
t =  (pi/2)*(x4 - L1-1-L2);
rt = (yc-y4);
xt = xc + rt.*cos(t);
yt = yc + L2 + rt.*sin(t);
mesh.p(1,ind4) = xt;
mesh.p(2,ind4) = yt;

t = (pi/2)*(X4 - L1-1-L2);
rt = (yc-Y4);
XT = X; XT(ine4) = xc + rt.*cos(t);
YT = Y; YT(ine4) = yc + L2 + rt.*sin(t);
mesh.dgnodes(:,1,:) = XT;
mesh.dgnodes(:,2,:) = YT;

% map x in [L1+1+L2+1 L1+1+L2+1+L3] -> x in [L1+R L1+R-L3] 
mesh.p(1,ind5) = L1 - (x5-(L1+1+L2+1));
% map y in [0 ymax] -> y in [ymax+R+L2+R+ymax ymax+R+L2+R] 
mesh.p(2,ind5) = ymax+R+L2+R+ymax - y5; 

XT = X; XT(ine5) = L1 - (X5-(L1+1+L2+1));
YT = Y; YT(ine5) = ymax+R+L2+R+ymax - Y5; 
mesh.dgnodes(:,1,:) = XT;
mesh.dgnodes(:,2,:) = YT;

mesh.p(1,:) = -mesh.p(1,:);
mesh.dgnodes(:,1,:) = - mesh.dgnodes(:,1,:);
% mesh.p = mesh.p/0.3048;
% mesh.dgnodes = mesh.dgnodes/0.3048;
plotboundary(mesh,3);
camlight headlight
lighting gouraud
material dull
camlight('right')
camlight('left')
lighting phong
set(gca,'AmbientLightColor',[0.6 0.6 0.6]*1.2)
set(gca,'FontSize',16);

