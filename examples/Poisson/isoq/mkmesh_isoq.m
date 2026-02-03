function mesh = mkmesh_isoq(porder)

[xl, xu] = isoq();

x1 = -0.04;
x2 = 0.013;

ind = (xu(:,1) <= x1);
xu1 = xu(ind,:);

ind = (xu(:,1) > x1);
xu2 = [xu1(end,:); xu(ind,:)];

ind = (xl(:,1) <= x2);
xl1 = xl(ind,:);

ind = (xl(:,1) > x2);
xl2 = [xl1(end,:); xl(ind,:)];

mesh1 = surfmesh2d(xl1, xu1, 64, 120, porder, [2.0 1.2], [5 0]);
mesh2 = surfmesh2d(xl2, xu2, 48, 120, porder, [2.0 1.5], [5 0]);

[mesh1, mesh2] = rightleft2d(mesh1, mesh2);

figure(1); clf; meshplot(mesh1);
hold on; meshplot(mesh2);
axis equal; axis tight;

mesh = mesh1;
[mesh.p, mesh.t] = connectmesh(mesh1.p', mesh1.t', mesh2.p', mesh2.t', 1e-5);
mesh.dgnodes = cat(3, mesh1.dgnodes, mesh2.dgnodes);
mesh.p = mesh.p';
mesh.t = mesh.t';

mesh.telem = mesh.tlocal;

figure(2); clf; meshplot(mesh);
axis on; axis equal; axis tight;
exportgraphics(gca,"mesh.png",'Resolution',200);

deltay = min(mesh.p(2,:));
L = max(mesh.p(1,:));
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(2,:)-deltay)<1e-6, @(p) p(1,:)> L-1e-4, @(p) ((p(1,:) < -1e-3) | (p(2,:) > 0.1)), @(p) abs(p(1,:))< 20 + 1e-6};
% axis symmetric, outflow, inflow, wall
mesh.boundarycondition = [6, 2, 1, 3]; % Set boundary condition for each boundary
mesh.f = facenumbering(mesh.p,mesh.t,1,mesh.boundaryexpr,[]);
mesh.periodicboundary = [];
mesh.periodicexpr = {};

figure(3); clf;
boundaryplot(mesh,1);
hold on;
for i = 2:4
  boundaryplot(mesh,i);
end
