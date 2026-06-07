function [mesh] = mkmesh_isoq2d2(porder, dR)

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

n1 = 48; m1 = 80; n2 = 36; m2 = 80;
mesh1 = surfmesh2d(xl1, xu1, n1, m1, porder, [2.0 1.2], [5 0]);
mesh2 = surfmesh2d(xl2, xu2, n2, m2, porder, [2.0 1.5], [5 0]);

[mesh1, mesh2] = rightleft2d(mesh1, mesh2);

figure(1); clf; meshplot(mesh1);
hold on; meshplot(mesh2);
axis equal; axis tight;

mesh = mesh1;
[mesh.p, mesh.t] = connectmesh(mesh1.p', mesh1.t', mesh2.p', mesh2.t', 1e-5);
mesh.dgnodes = cat(3, mesh1.dgnodes, mesh2.dgnodes);
mesh.p = mesh.p';
mesh.t = mesh.t';

mesh.p(2,:) = mesh.p(2,:) + dR;
mesh.dgnodes(:,2,:) = mesh.dgnodes(:,2,:) + dR;

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
mesh.dist = meshdist3(mesh.f,mesh.dgnodes,mesh.perm,4); % distance to the wall

figure(3); clf;
boundaryplot(mesh,1);
hold on;
for i = 2:4
  boundaryplot(mesh,i);
end
