function [mesh, xc, xd] = mkmesh_isoq2d(porder)

[xl, xu] = isoq();

ind = xl(:,2) >= 0.012;
xl = xl(ind,:);

ind = xu(:,2) >= 0.015;
xu = xu(ind,:);

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

n1 = 48/2; m1 = 120/8;
n2 = 48/2; m2 = 120/8;
mesh1 = surfmesh2d(xl1, xu1, n1, m1, porder, [1.0 1.2], [2 0]);
mesh2 = surfmesh2d(xl2, xu2, n2, m2, porder, [2.0 1.5], [2 0]);

xc = mesh1.p(:,1:n1+1:end);
x = mesh1.dgnodes(:,1,:); x=x(:);
y = mesh1.dgnodes(:,2,:); y=y(:);
a = xc(:,1); b = xc(:,end);
ind = y < a(2) + (b(2)-a(2))*(x-a(1))/(b(1)-a(1)) + 1e-6;
xd = [x(ind) y(ind)]';

[mesh1, mesh2] = rightleft2d(mesh1, mesh2);

figure(1); clf; meshplot(mesh1);
hold on; meshplot(mesh2);
axis equal; % axis tight;

mesh = mesh1;
[mesh.p, mesh.t] = connectmesh(mesh1.p', mesh1.t', mesh2.p', mesh2.t', 1e-5);
mesh.dgnodes = cat(3, mesh1.dgnodes, mesh2.dgnodes);
mesh.p = mesh.p';
mesh.t = mesh.t';
mesh.telem = mesh.tlocal;

figure(2); clf; meshplot(mesh);
axis on; axis equal; %axis tight;
%exportgraphics(gca,"mesh.png",'Resolution',200);
hold on;
plot(xc(1,:),xc(2,:),'o');
plot(xd(1,:),xd(2,:),'*');

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
plot(xc(1,:),xc(2,:),'o');
plot(xd(1,:),xd(2,:),'*');