function mesh = mkmesh_orion3(porder, D)

if nargin<2
  D = 0.12;
end

load shockcurve.mat;
[xl, xs] = orion(D);

ymin = 0.0;
yrear = 0.5;
yshock = 7.2;
xmaxxl = max(xl(:,1));
x1 = -0.25;
x2 = 1.05;
x3 = 0.4812;
x4 = 0.8477;
x5 = 5.0;
x6 = 9;

t = linspace(pi, pi/2, 2000)';
xc = 13 + 14.4*cos(t);
yc = ymin + 8.8*sin(t);
ind = xc <=x6;
xu = [xc(ind) yc(ind)];

ind = (xm(:,1) <= x6);
xm = xm(ind,:);

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

xm1 = [0 0];
j = 1;
for i = 1:size(xm,1)
  [~, ~, position] = find_perpendicular_point(xl1(end,1), xl1(end,2), xu1(end,1), xu1(end,2), xm(i,1), xm(i,2));
  if position<=0
    xm1(j,:) = xm(i,:);
    j = j + 1;
  end
end

xm2 = [0 0];
j = 1;
for i = 1:size(xm,1)
  [~, ~, position1] = find_perpendicular_point(xl1(end,1), xl1(end,2), xu1(end,1), xu1(end,2), xm(i,1), xm(i,2));  
  [~, ~, position2] = find_perpendicular_point(xl2(end,1), xl2(end,2), xu2(end,1), xu2(end,2), xm(i,1), xm(i,2));
  if (position1>=0) && (position2>=0)
    xm2(j,:) = xm(i,:);
    j = j + 1;
  end
end

xm3 = [0 0];
j = 1;
for i = 1:size(xm,1)
  [~, ~, position] = find_perpendicular_point(xl2(end,1), xl2(end,2), xu2(end,1), xu2(end,2), xm(i,1), xm(i,2));
  if position<=0
    xm3(j,:) = xm(i,:);
    j = j + 1;
  end
end

mesh1 = surfmesh2d_shockcurve(xl1, xu1, xm1, 40, 80, 15, porder, [2.5 2], [8 7], [5 0], 0);
mesh2 = surfmesh2d_shockcurve(xl2, xu2, xm2, 16, 80, 15, porder, [.1 .1], [8 7], [5 0], 0);
mesh3 = surfmesh2d_shockcurve(xl3, xu3, xm3, 30, 80, 15, porder, [3 2], [8 7], [5 0], 0.65);

mesh3.dgnodes = fixdgnodes(mesh3.p,mesh3.dgnodes,1e-5);

[mesh1, mesh2] = rightleft2d(mesh1, mesh2);
[mesh2, mesh3] = rightleft2d(mesh2, mesh3);

xi = linspace(0,1,200)';
XA = mesh3.cgright(64,:);
XB = [x6 yshock];
x = XA(1,1) + (XB(1,1)-XA(1,1))*xi;
y = XA(1,2) + (XB(1,2)-XA(1,2))*xi;
xl4 = [x y];

XA = mesh3.cgright(end,:);
ind = (xu(:,1) > XA(1));
xu4 = [XA; xu(ind,:)];
xu4(end,1) = x6;

cgw = mesh3.cgright;
dgw = mesh3.dgright;
ind = (cgw(:,1) > xl4(1,1));
cga = [xl4(1,:); cgw(ind,:)];
ind = (dgw(:,1) > xl4(1,1));
dga = [xl4(1,:); dgw(ind,:)];

mesh4 = surfmesh2d(xl4, xu4, 30, cga, porder, [0 0], [0 0], 0.25, [], dga);

ind = (cgw(:,1) < xl4(1,1));
cgb = [cgw(ind,:); xl4(1,:)];

xi = logdec(linspace(0,1,161)',4.2);
XA = [x6 xl3(end,2)];
XB = xl3(end,:);
x = XA(1,1) + (XB(1,1)-XA(1,1))*xi;
y = XA(1,2) + (XB(1,2)-XA(1,2))*xi;
cgc = [x y];

xi = loginc(logdec(linspace(0,1,30)',4),3);
XA = [x6 yshock];
XB = [x6 xl3(end,2)];
x = XA(1,1) + (XB(1,1)-XA(1,1))*xi;
y = XA(1,2) + (XB(1,2)-XA(1,2))*xi;
cgd = [x y];

cga = mesh4.cgbottom;
cg5 = [cgd; cgc(2:end,:); cgb(2:end,:); cga(2:end-1,:)];
cg5 = cg5(end:-1:1,:);
poly2gmsh('domain5.geo', cg5, 0.6);
gmshmatlab('domain5', '-2 -format msh3');
[p1,t1] = gmsh2ptnew('domain5.msh',2,1);
p1 = p1'; t1 = t1';
mesh5 = mkmesh(p1,t1,porder,{'true'}, 1,1);
mesh5.p = p1';
mesh5.t = t1';

ind = abs(p1(:,2) - xl3(end,2)) <= 1e-4;
cga = p1(ind,:);
[~,ind] = sort(cga(:,1));
cga = cga(ind,:);

xi = logdec(linspace(0,1,6)',0.5);
XA = [x6 xl3(end,2)];
XB = [x6 yrear];
x = XA(1,1) + (XB(1,1)-XA(1,1))*xi;
y = XA(1,2) + (XB(1,2)-XA(1,2))*xi;
cgb = [x y];

xi = linspace(0,1,51)';
XA = [x6 yrear];
XB = [xl3(end,1) yrear];
x = XA(1,1) + (XB(1,1)-XA(1,1))*xi;
y = XA(1,2) + (XB(1,2)-XA(1,2))*xi;
cgc = [x y];

xi = logdec(linspace(0,1,20)',3.6);
XA = [xl3(end,1) yrear];
XB = xl3(end,:);
x = XA(1,1) + (XB(1,1)-XA(1,1))*xi;
y = XA(1,2) + (XB(1,2)-XA(1,2))*xi;
cgd = [x y];

cg6 = [cga; cgb(2:end,:); cgc(2:end,:); cgd(2:end-1,:)];
cg6 = cg6(end:-1:1,:);
poly2gmsh('domain6.geo', cg6, 0.3);
gmshmatlab('domain6', '-2 -format msh3');
[p1,t1] = gmsh2ptnew('domain6.msh',2,1);
p1 = p1'; t1 = t1';
mesh6 = mkmesh(p1,t1,porder,{'true'}, 1,1);
mesh6.p = p1';
mesh6.t = t1';

ind = abs(p1(:,2) - yrear) <= 1e-4;
cga = p1(ind,:);
[~,ind] = sort(cga(:,1));
cga = cga(ind,:);
cgb = linspace(ymin, yrear, 7)';

mesh7 = mkmesh_cartesian(cga(:,1),cgb,porder,1,1,1);

XA = mesh3.cgright(1,:);
XB = mesh3.cgright(end,:);
expr = @(p) abs((p(1)-XA(1))/(XB(1)-XA(1))-(p(2)-XA(2))/(XB(2)-XA(2)))<1e-5;
[mesh3.p, mesh5.p] = fixp(mesh3.p, mesh5.p, 2, 1e-4, expr);
[mesh3.dgnodes, mesh5.dgnodes] = fixp(mesh3.dgnodes, mesh5.dgnodes, 2, 1e-2, expr);

XA = xl4(1,:);
XB = xl4(end,:);
expr = @(p) abs((p(1)-XA(1))/(XB(1)-XA(1))-(p(2)-XA(2))/(XB(2)-XA(2)))<1e-5;
[mesh4.p, mesh5.p] = fixp(mesh4.p, mesh5.p, 2, 1e-4, expr);
[mesh4.dgnodes, mesh5.dgnodes] = fixp(mesh4.dgnodes, mesh5.dgnodes, 2, 1e-3, expr);

mesh5.xpe = mesh1.xpe;
mesh6.xpe = mesh1.xpe;
mesh7.xpe = mesh1.xpe;

figure(1); clf; 
hold on;
meshplot(mesh1);
meshplot(mesh2);
meshplot(mesh3);
meshplot(mesh4);
meshplot(mesh5);
meshplot(mesh6);
meshplot(mesh7);

% x = mesh1.dgnodes(:,1,:); y = mesh1.dgnodes(:,2,:); plot(x(:), y(:), 'ob');
% x = mesh2.dgnodes(:,1,:); y = mesh2.dgnodes(:,2,:); plot(x(:), y(:), 'or');
% x = mesh3.dgnodes(:,1,:); y = mesh3.dgnodes(:,2,:); plot(x(:), y(:), 'ob');
% x = mesh4.dgnodes(:,1,:); y = mesh4.dgnodes(:,2,:); plot(x(:), y(:), 'ok');
% x = mesh5.dgnodes(:,1,:); y = mesh5.dgnodes(:,2,:); plot(x(:), y(:), 'or');
% x = mesh6.dgnodes(:,1,:); y = mesh6.dgnodes(:,2,:); plot(x(:), y(:), 'ob');
% x = mesh7.dgnodes(:,1,:); y = mesh7.dgnodes(:,2,:); plot(x(:), y(:), 'or');
% x = mesh3.p(1,:); y = mesh3.p(2,:); plot(x(:), y(:), 'ob', 'MarkerSize', 6, 'LineWidth', 2);
% x = mesh5.p(1,:); y = mesh5.p(2,:); plot(x(:), y(:), 'or', 'MarkerSize', 6, 'LineWidth', 2);

mesh = mesh1;
[mesh.p, mesh.t] = connectmesh(mesh1.p', mesh1.t', mesh2.p', mesh2.t', 1e-5);
[mesh.p, mesh.t] = connectmesh(mesh.p, mesh.t, mesh3.p', mesh3.t', 1e-5);
[mesh.p, mesh.t] = connectmesh(mesh.p, mesh.t, mesh4.p', mesh4.t', 1e-5);
[mesh.p, mesh.t] = connectmesh(mesh.p, mesh.t, mesh5.p', mesh5.t', 1e-5);
[mesh.p, mesh.t] = connectmesh(mesh.p, mesh.t, mesh6.p', mesh6.t', 1e-5);
[mesh.p, mesh.t] = connectmesh(mesh.p, mesh.t, mesh7.p', mesh7.t', 1e-5);
mesh.dgnodes = cat(3, mesh1.dgnodes, mesh2.dgnodes);
mesh.dgnodes = cat(3, mesh.dgnodes, mesh3.dgnodes);
mesh.dgnodes = cat(3, mesh.dgnodes, mesh4.dgnodes);
mesh.dgnodes = cat(3, mesh.dgnodes, mesh5.dgnodes);
mesh.dgnodes = cat(3, mesh.dgnodes, mesh6.dgnodes);
mesh.dgnodes = cat(3, mesh.dgnodes, mesh7.dgnodes);

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

figure(2); clf; 
meshplot(mesh);


% xi = linspace(0,1,200)';
% XA = xl3(end,:);
% XB = [xl3(end,1) ymin];
% x = XA(1,1) + (XB(1,1)-XA(1,1))*xi;
% y = XA(1,2) + (XB(1,2)-XA(1,2))*xi;
% xl4 = [x y];
% 
% XA = [x6 yrear];
% XB = [x6 ymin];
% x = XA(1,1) + (XB(1,1)-XA(1,1))*xi;
% y = XA(1,2) + (XB(1,2)-XA(1,2))*xi;
% xu4 = [x y];
% 
% mesh4 = surfmesh2d(xl4, xu4, 20, 30, porder, [2 0], [3 3]);
% figure(1); clf; meshplot(mesh4);
% 
% xl5 = mesh4.cgleft;
% dg5 = mesh4.dgleft;
% 
% XA = mesh3.cgright(end,:);
% ind = (xu(:,1) > XA(1));
% xu5 = [XA; xu(ind,:)];
% xu5(end,1) = x6;
% 
% xw = mesh3.cgright;
% dgw = mesh3.dgright;
% 
% hold on; 
% plot(xl5(:,1), xl5(:,2), '-o');
% plot(xu5(:,1), xu5(:,2), '-o');
% 
% mesh5 = surfmesh2d(xl5, xu5, xl5, xw, porder, [0 0], [0 0], 0.5, dg5, dgw);
% meshplot(mesh5);
% 
% mesh = mesh1;
% [mesh.p, mesh.t] = connectmesh(mesh1.p', mesh1.t', mesh2.p', mesh2.t', 1e-5);
% [mesh.p, mesh.t] = connectmesh(mesh.p, mesh.t, mesh3.p', mesh3.t', 1e-5);
% [mesh.p, mesh.t] = connectmesh(mesh.p, mesh.t, mesh4.p', mesh4.t', 1e-5);
% [mesh.p, mesh.t] = connectmesh(mesh.p, mesh.t, mesh5.p', mesh5.t', 1e-5);
% mesh.dgnodes = cat(3, mesh1.dgnodes, mesh2.dgnodes);
% mesh.dgnodes = cat(3, mesh.dgnodes, mesh3.dgnodes);
% mesh.dgnodes = cat(3, mesh.dgnodes, mesh4.dgnodes);
% mesh.dgnodes = cat(3, mesh.dgnodes, mesh5.dgnodes);
% 
% mesh.p = mesh.p';
% mesh.t = mesh.t';
% deltay = min(mesh.p(2,:));
% L = max(mesh.p(1,:));
% % expressions for domain boundaries
% mesh.boundaryexpr = {@(p) abs(p(2,:)-deltay)<1e-6, @(p) p(1,:)> L-1e-2, @(p) ((p(1,:) < -1e-3) | (p(2,:) > 2.6)), @(p) abs(p(1,:))< 20 + 1e-6};
% % axis symmetric, outflow, inflow, wall
% mesh.boundarycondition = [6, 2, 1, 3]; % Set boundary condition for each boundary
% mesh.f = facenumbering(mesh.p,mesh.t,1,mesh.boundaryexpr,[]);
% mesh.periodicboundary = [];
% mesh.periodicexpr = {};
% 
% 
% figure(2); clf; 
% meshplot(mesh);
% hold on;
% % meshplot(meshht);
% % size(xs)
% % size(xf)
% % plot(meshht.cgmiddle(:,1), meshht.cgmiddle(:,2), '-r', 'LineWidth',1);
% % x = mesh1.dgnodes(:,1,:); y = mesh1.dgnodes(:,2,:); plot(x(:), y(:), 'ob');
% % x = mesh2.dgnodes(:,1,:); y = mesh2.dgnodes(:,2,:); plot(x(:), y(:), 'or');
% % x = mesh3.dgnodes(:,1,:); y = mesh3.dgnodes(:,2,:); plot(x(:), y(:), 'ob');
% % x = mesh4.dgnodes(:,1,:); y = mesh4.dgnodes(:,2,:); plot(x(:), y(:), 'or');
% % x = mesh5.dgnodes(:,1,:); y = mesh5.dgnodes(:,2,:); plot(x(:), y(:), 'ok');
% plot(mesh1.dgbottom(:,1), mesh1.dgbottom(:,2), '-r', 'LineWidth', 1.0);
% % plot(mesh1.dgtop(:,1), mesh1.dgtop(:,2), '-r',  'LineWidth', 1.0);
% % plot(mesh1.dgleft(:,1), mesh1.dgleft(:,2), '-r',  'LineWidth', 1.0);
% % plot(mesh1.dgright(:,1), mesh1.dgright(:,2), '-r',  'LineWidth', 1.0);
% plot(mesh2.dgbottom(:,1), mesh2.dgbottom(:,2), '-r', 'LineWidth', 1.0);
% % plot(mesh2.dgtop(:,1), mesh2.dgtop(:,2), '-r',  'LineWidth', 1.0);
% % plot(mesh2.dgleft(:,1), mesh2.dgleft(:,2), '-r',  'LineWidth', 1.0);
% % plot(mesh2.dgright(:,1), mesh2.dgright(:,2), '-r',  'LineWidth', 1.0);
% plot(mesh3.dgbottom(:,1), mesh3.dgbottom(:,2), '-r', 'LineWidth', 1.0);
% % plot(mesh3.dgtop(:,1), mesh3.dgtop(:,2), '-r',  'LineWidth', 1.0);
% % plot(mesh3.dgleft(:,1), mesh3.dgleft(:,2), '-r',  'LineWidth', 1.0);
% % plot(mesh3.dgright(:,1), mesh3.dgright(:,2), '-r',  'LineWidth', 1.0);
% plot(mesh4.dgbottom(:,1), mesh4.dgbottom(:,2), '-r', 'LineWidth', 1.0);
% % plot(mesh4.dgtop(:,1), mesh4.dgtop(:,2), '-r',  'LineWidth', 1.0);
% % plot(mesh4.dgleft(:,1), mesh4.dgleft(:,2), '-r',  'LineWidth', 1.0);
% % plot(mesh4.dgright(:,1), mesh4.dgright(:,2), '-r',  'LineWidth', 1.0);
% % plot(mesh5.dgbottom(:,1), mesh5.dgbottom(:,2), '-r', 'LineWidth', 1.0);
% % plot(mesh5.dgtop(:,1), mesh5.dgtop(:,2), '-r',  'LineWidth', 1.0);
% % plot(mesh5.dgleft(:,1), mesh5.dgleft(:,2), '-r',  'LineWidth', 1.0);
% % plot(mesh5.dgright(:,1), mesh5.dgright(:,2), '-r',  'LineWidth', 1.0);
% set(gca,'FontSize',20); 
% set(gca,'LooseInset',get(gca,'TightInset'))
% axis equal; axis on; box on;
% xlabel("$z$", 'interpreter', 'latex', 'FontSize', 28);
% ylabel("$r$", 'interpreter', 'latex', 'FontSize', 28);
% 
% 
