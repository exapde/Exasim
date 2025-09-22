function [mesh, meshht] = mkmesh_orion(porder, D)

if nargin<2
  D = 0.12;
end

load shockcurve.mat;
[xl, xs] = orion(D);

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

mesh1 = surfmesh2d_shockcurve(xl1, xu1, xm1, 48, 80, 15, porder, [2.5 2], [8 7], [5 0], 0);
mesh2 = surfmesh2d_shockcurve(xl2, xu2, xm2, 20, 80, 15, porder, [.1 .1], [8 7], [5 0], 0);
mesh3 = surfmesh2d_shockcurve(xl3, xu3, xm3, 24, 80, 15, porder, [3 1], [8 7], [5 0], 0.65);

[mesh1, mesh2] = rightleft2d(mesh1, mesh2);
[mesh2, mesh3] = rightleft2d(mesh2, mesh3);

% max(abs([mesh1.cgright-mesh2.cgleft]))
% max(abs([mesh1.dgright-mesh2.dgleft]))
% max(abs([mesh2.cgright-mesh3.cgleft]))
% max(abs([mesh2.dgright-mesh3.dgleft]))
  
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

mesh4 = surfmesh2d(xl4, xu4, 20, xr, porder, [0 0], [0 0], 0.5, [], dgr);

xl5 = mesh4.cgtop;
dg5 = mesh4.dgtop;

ind = xw(:,2) > 1.1;
xq = [xr(end,:); xw(ind,:)];
dgq = reshape(dgw(:,(nr+1):end,:), [], 2);

XA = xq(end,:);
% XB = [7 xq(end,2)+0.4];
% x = XA(1,1) + (XB(1,1)-XA(1,1))*xi;
% y = XA(1,2) + (XB(1,2)-XA(1,2))*xi;
% xu5 = [x y];
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

n = 12;
if (D==0.2)
  xb = [0.665 0.745 max(xs(:,1))];
elseif (D==0.15)
  xb = [0.625 0.77 max(xs(:,1))];
elseif (D==0.12)
  xb = [0.59 0.79 max(xs(:,1))];    
elseif (D==0.1)
  xb = [0.57 0.8 max(xs(:,1))];  
else
  error("D can only be 0.1, 0.12, 0.15, or 0.2");
end

xl6 = mesh1.cgbottom;
dg6 = mesh1.dgbottom;
ind = xs(:,1) <= xb(1);
xu6 = xs(ind,:);
mesh6 = surfmesh2d(xl6, xu6, xl6, n, porder, [0 0], [4 0], 0, dg6, []);

xl7 = mesh2.cgbottom;
dg7 = mesh2.dgbottom;
ind = (xb(1) <= xs(:,1)) & (xs(:,1) <= xb(2));

xu7 = [xu6(end,:); xs(ind,:)];
mesh7 = surfmesh2d(xl7, xu7, xl7, n, porder, [0 0], [4 0], 0, dg7, []);

xl8 = mesh3.cgbottom;
dg8 = mesh3.dgbottom;
ind = abs(xs(:,1) - xb(3)) < 1e-8;
xu9 = xs(ind,:);
ind = (xb(2) <= xs(:,1)) & (xs(:,1) < xb(3));
xu8 = [xu7(end,:); xs(ind,:); xu9(1,:)];
mesh8 = surfmesh2d(xl8, xu8, xl8, n, porder, [0 0], [4 0], 0, dg8, []);

xl9 = mesh4.cgbottom;
dg9 = mesh4.dgbottom;
mesh9 = surfmesh2d(xl9, xu9, xl9, n, porder, [0 0], [4 0], 0, dg9, []);

meshht = mesh;
[meshht.p, meshht.t] = connectmesh(mesh6.p', mesh6.t', mesh7.p', mesh7.t', 1e-5);
[meshht.p, meshht.t] = connectmesh(meshht.p, meshht.t, mesh8.p', mesh8.t', 1e-5);
[meshht.p, meshht.t] = connectmesh(meshht.p, meshht.t, mesh9.p', mesh9.t', 1e-5);
meshht.dgnodes = cat(3, mesh6.dgnodes, mesh7.dgnodes);
meshht.dgnodes = cat(3, meshht.dgnodes, mesh8.dgnodes);
meshht.dgnodes = cat(3, meshht.dgnodes, mesh9.dgnodes);
meshht.p = meshht.p';
meshht.t = meshht.t';

meshht.t = meshht.t([1 4 3 2],:);
ind = reshape(1:((porder+1)*(porder+1)), [(porder+1) (porder+1)])';
meshht.dgnodes = meshht.dgnodes(ind(:),:,:);

p6 = 0.5*(mesh6.cgbottom + mesh6.cgtop);
p7 = 0.5*(mesh7.cgbottom + mesh7.cgtop);
p8 = 0.5*(mesh8.cgbottom + mesh8.cgtop);
p9 = 0.5*(mesh9.cgbottom + mesh9.cgtop);
mid = [p6; p7(2:end,:); p8(2:end,:); p9(2:end,:)];
save midcurve.mat mid;

meshht.boundaryexpr = {@(p) abs(p(2,:)-deltay)<1e-4, @outerwall, @(p) abs(p(1,:))< 20 + 1e-6};
meshht.boundarycondition = [1, 2, 3]; % Set boundary condition for each boundary
meshht.f = facenumbering(meshht.p,meshht.t,1,meshht.boundaryexpr,[]);
meshht.periodicboundary = [];
meshht.periodicexpr = {};

figure(1);clf;boundaryplot(meshht,1);
hold on; boundaryplot(meshht,2);
boundaryplot(meshht,3);
plot(mid(:,1), mid(:,2), 'ok');

figure(2); clf; 
meshplot(mesh);
hold on;
meshplot(meshht);
set(gca,'FontSize',20); 
set(gca,'LooseInset',get(gca,'TightInset'))
axis equal; axis tight; axis on; box on;
xlabel("$z$", 'interpreter', 'latex', 'FontSize', 28);
ylabel("$r$", 'interpreter', 'latex', 'FontSize', 28);

figure(3); clf; 
meshplot(mesh);
hold on;
meshplot(meshht);
% size(xs)
% size(xf)
% plot(meshht.cgmiddle(:,1), meshht.cgmiddle(:,2), '-r', 'LineWidth',1);
% x = mesh1.dgnodes(:,1,:); y = mesh1.dgnodes(:,2,:); plot(x(:), y(:), 'ob');
% x = mesh2.dgnodes(:,1,:); y = mesh2.dgnodes(:,2,:); plot(x(:), y(:), 'or');
% x = mesh3.dgnodes(:,1,:); y = mesh3.dgnodes(:,2,:); plot(x(:), y(:), 'ob');
% x = mesh4.dgnodes(:,1,:); y = mesh4.dgnodes(:,2,:); plot(x(:), y(:), 'or');
% x = mesh5.dgnodes(:,1,:); y = mesh5.dgnodes(:,2,:); plot(x(:), y(:), 'ok');
plot(mesh1.dgbottom(:,1), mesh1.dgbottom(:,2), '-r', 'LineWidth', 1.0);
% plot(mesh1.dgtop(:,1), mesh1.dgtop(:,2), '-r',  'LineWidth', 1.0);
% plot(mesh1.dgleft(:,1), mesh1.dgleft(:,2), '-r',  'LineWidth', 1.0);
% plot(mesh1.dgright(:,1), mesh1.dgright(:,2), '-r',  'LineWidth', 1.0);
plot(mesh2.dgbottom(:,1), mesh2.dgbottom(:,2), '-r', 'LineWidth', 1.0);
% plot(mesh2.dgtop(:,1), mesh2.dgtop(:,2), '-r',  'LineWidth', 1.0);
% plot(mesh2.dgleft(:,1), mesh2.dgleft(:,2), '-r',  'LineWidth', 1.0);
% plot(mesh2.dgright(:,1), mesh2.dgright(:,2), '-r',  'LineWidth', 1.0);
plot(mesh3.dgbottom(:,1), mesh3.dgbottom(:,2), '-r', 'LineWidth', 1.0);
% plot(mesh3.dgtop(:,1), mesh3.dgtop(:,2), '-r',  'LineWidth', 1.0);
% plot(mesh3.dgleft(:,1), mesh3.dgleft(:,2), '-r',  'LineWidth', 1.0);
% plot(mesh3.dgright(:,1), mesh3.dgright(:,2), '-r',  'LineWidth', 1.0);
plot(mesh4.dgbottom(:,1), mesh4.dgbottom(:,2), '-r', 'LineWidth', 1.0);
% plot(mesh4.dgtop(:,1), mesh4.dgtop(:,2), '-r',  'LineWidth', 1.0);
% plot(mesh4.dgleft(:,1), mesh4.dgleft(:,2), '-r',  'LineWidth', 1.0);
% plot(mesh4.dgright(:,1), mesh4.dgright(:,2), '-r',  'LineWidth', 1.0);
% plot(mesh5.dgbottom(:,1), mesh5.dgbottom(:,2), '-r', 'LineWidth', 1.0);
% plot(mesh5.dgtop(:,1), mesh5.dgtop(:,2), '-r',  'LineWidth', 1.0);
% plot(mesh5.dgleft(:,1), mesh5.dgleft(:,2), '-r',  'LineWidth', 1.0);
% plot(mesh5.dgright(:,1), mesh5.dgright(:,2), '-r',  'LineWidth', 1.0);
set(gca,'FontSize',20); 
set(gca,'LooseInset',get(gca,'TightInset'))
axis equal; axis on; box on;
axis([-0.2 xb(end)+3*D 0 max(xl(:,2))+2*D]);
xlabel("$z$", 'interpreter', 'latex', 'FontSize', 28);
ylabel("$r$", 'interpreter', 'latex', 'FontSize', 28);

