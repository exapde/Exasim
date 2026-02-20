function mesh = mkmesh_sharpb2(porder)

[xl, xu, T1] = sharpb2smooth();

x1 = T1(1);
x2 = 0.0;

ind = (xl(:,1) <= x1);
xl1 = xl(ind,:);

ind = (xu(:,1) <= x2);
xu1 = xu(ind,:);

ind = (xl(:,1) > x1) & (xl(:,1) <= 1.25);
xl2 = [xl1(end,:); xl(ind,:)];

ind = (xu(:,1) > x2) & (xu(:,1) <= 1.22);
xu2 = [xu1(end,:); xu(ind,:)];

mesh1 = surfmesh2d(xl1, xu1, 36, 100, porder, [2.0 1.2], [5 0]);
mesh2 = surfmesh2d(xl2, xu2, 120, 100, porder, [3.5 1.5], [5 0]);

[mesh1, mesh2] = rightleft2d(mesh1, mesh2);

figure(1); clf; meshplot(mesh1);
axis equal; axis tight;

figure(2); clf; meshplot(mesh2);
axis equal; axis tight;

mesh = mesh1;
[mesh.p, mesh.t] = connectmesh(mesh1.p', mesh1.t', mesh2.p', mesh2.t', 1e-5);
mesh.dgnodes = cat(3, mesh1.dgnodes, mesh2.dgnodes);
mesh.p = mesh.p';
mesh.t = mesh.t';

mesh.p(2,:) = mesh.p(2,:) + 1e-4;
mesh.dgnodes(:,2,:) = mesh.dgnodes(:,2,:) + 1e-4;

mesh.telem = mesh.tlocal;
figure(3); clf; meshplot(mesh);
axis on; axis equal; axis tight;
exportgraphics(gca,"mesh.png",'Resolution',200);

xmin = min(mesh.p(1,:));
ymin = min(mesh.p(2,:));
[y2,ind] = max(mesh.p(2,:));
x2 = mesh.p(1,ind);
[x3,ind] = max(mesh.p(1,:));
y3 = mesh.p(2,ind);
x4 = 0;
y4 = 0.04;
% expressions for domain boundaries
%mesh.boundaryexpr = {@(p) abs(p(2,:)-ymin)<1e-6, @(p) p(2,:) > ymin + (y2-ymin)/(x2-xmin)*(p(1,:)-xmin) - 1e-4, @(p) p(2,:) < y3 + 1e-4, @(p) abs(p(1,:))< 20 + 1e-6};
mesh.boundaryexpr = {@(p) abs(p(2,:)-ymin)<1e-6, @(p) p(2,:) > ymin + (y4-ymin)/(x4-xmin)*(p(1,:)-xmin) - 1e-4,  @(p) p(2,:) > y4 + (y2-y4)/(x2-x4)*(p(1,:)-x4) - 1e-4, @(p) p(2,:) < y3 + 1e-4, @(p) abs(p(1,:))< 20 + 1e-6};
% axis symmetric, outflow, inflow, wall
mesh.boundarycondition = [6, 2, 1, 3]; % Set boundary condition for each boundary
mesh.f = facenumbering(mesh.p,mesh.t,1,mesh.boundaryexpr,[]);
mesh.periodicboundary = [];
mesh.periodicexpr = {};

figure(4); clf;
boundaryplot(mesh,1,'b');
hold on;
colors = ['r', 'g', 'k', 'y'];
for i = 2:5
  boundaryplot(mesh,i, colors(i-1));
end
axis equal;
%axis tight;
