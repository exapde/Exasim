function mesh = mkmesh_orion6(porder)

xl = orion();

ymin = 0.0;
xmaxxl = max(xl(:,1));
x1 = -0.25;
x2 = 1.05;
x3 = 0.4812;
x4 = 0.8477;
x5 = 5.6;
x6 = 9;

t = linspace(pi, pi/2, 2000)';
xc = 13 + 14.4*cos(t);
yc = ymin + 8.8*sin(t);
ind = xc <=x6;
xu = [xc(ind) yc(ind)];

ind = (xl(:,1) <= xmaxxl);
xl = xl(ind,:);

ind = (xu(:,1) <= x6);
xu = xu(ind,:);

%figure(1);clf;plot(xl(:,1),xl(:,2),xu(:,1),xu(:,2),xm(:,1),xm(:,2));

ind = (xu(:,1) <= x1);
xu1 = xu(ind,:);

ind = (xu(:,1) > x1 & xu(:,1) <= x2);
xu2 = [xu1(end,:); xu(ind,:)];

ind = (xu(:,1) > x2 & xu(:,1) <= x5);
xu3 = [xu2(end,:); xu(ind,:)];

ind = (xl(:,1) <= x3);
xl1 = xl(ind,:);

ind = (xl(:,1) > x3 & xl(:,1) <= x4);
xl2 = [xl1(end,:); xl(ind,:)];

ind = (xl(:,1) > x4);
xl3 = [xl2(end,:); xl(ind,:)];

mesh1 = surfmesh2d(xl1, xu1, 40, 96, porder, [2.5 2], [7 0]);
mesh2 = surfmesh2d(xl2, xu2, 16, 96, porder, [0.2 0.1], [7 0]);
mesh3 = surfmesh2d(xl3, xu3, 28, 96, porder, [3 2], [7 0], 0.5);

[mesh1, mesh2] = rightleft2d(mesh1, mesh2);
[mesh2, mesh3] = rightleft2d(mesh2, mesh3);

xw = mesh3.cgright;
dgw = mesh3.dgright;

ind = xw(:,2) <= 1.1;
xr = xw(ind,:);
nr = size(xr,1) - 1;

ne  = size(xw,1) - 1;
dgw = reshape(dgw, porder+1, ne, 2);
dgr = reshape(dgw(:,1:nr,:), [(porder+1)*nr 2]);

xi = linspace(0,1,200)';
XA = xl3(end,:);
XB = [xl3(end,1) ymin];
x = XA(1,1) + (XB(1,1)-XA(1,1))*xi;
y = XA(1,2) + (XB(1,2)-XA(1,2))*xi;
xl4 = [x y];

XA = xr(end,:);
XB = [x6 ymin];
x = XA(1,1) + (XB(1,1)-XA(1,1))*xi;
y = XA(1,2) + (XB(1,2)-XA(1,2))*xi;
xu4 = [x y];

mesh4 = surfmesh2d(xl4, xu4, 24, xr, porder, [0 0], [0 0], 0.5, [], dgr);

xl5 = mesh4.cgtop;
dg5 = mesh4.dgtop;

ind = xw(:,2) > 1.1;
xq = [xr(end,:); xw(ind,:)];
dgq = reshape(dgw(:,(nr+1):end,:), [], 2);

XA = xq(end,:);
ind = (xu(:,1) > XA(1));
xu5 = [XA; xu(ind,:)];
xu5(end,1) = x6;
mesh5 = surfmesh2d(xl5, xu5, xl5, xq, porder, [0 0], [0 0], 0.5, dg5, dgq);

mesh = mesh1;
[mesh.p, mesh.t] = connectmesh(mesh1.p', mesh1.t', mesh2.p', mesh2.t', 1e-5);
[mesh.p, mesh.t] = connectmesh(mesh.p, mesh.t, mesh3.p', mesh3.t', 1e-5);
[mesh.p, mesh.t] = connectmesh(mesh.p, mesh.t, mesh4.p', mesh4.t', 1e-5);
[mesh.p, mesh.t] = connectmesh(mesh.p, mesh.t, mesh5.p', mesh5.t', 1e-5);
mesh.dgnodes = cat(3, mesh1.dgnodes, mesh2.dgnodes);
mesh.dgnodes = cat(3, mesh.dgnodes, mesh3.dgnodes);
mesh.dgnodes = cat(3, mesh.dgnodes, mesh4.dgnodes);
mesh.dgnodes = cat(3, mesh.dgnodes, mesh5.dgnodes);

mesh.p = mesh.p';
mesh.t = mesh.t';
deltay = min(mesh.p(2,:));
L = max(mesh.p(1,:));
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(2,:)-deltay)<1e-6, @(p) p(1,:)> L-1e-2, @(p) ((p(1,:) < -1e-3) | (p(2,:) > 2.6)), @(p) abs(p(1,:))< 20 + 1e-6};
% axis symmetric, outflow, inflow, wall
mesh.boundarycondition = [6, 2, 1, 3]; % Set boundary condition for each boundary
mesh.f = facenumbering(mesh.p,mesh.t,1,mesh.boundaryexpr,[]);
mesh.periodicboundary = [];
mesh.periodicexpr = {};

figure(1); clf;
hold on;
meshplot(mesh);

