function [mesh, rho, drhodx, drhody] = mkmesh_square(n, porder, elemtype)

[mesh.p,mesh.t] = squaremesh(n,n,1,elemtype);
% mesh.p(1,:) = loginc(mesh.p(1,:),2);
% mesh.p(2,:) = loginc(mesh.p(2,:),2);

mesh.p(1,:) = 2*mesh.p(1,:) - 1;
mesh.p(2,:) = 2*mesh.p(2,:) - 1;

% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(2,:)+1)<1e-8, @(p) abs(p(1,:)-1)<1e-8, @(p) abs(p(2,:)-1)<1e-8, @(p) abs(p(1,:)+1)<1e-8};
mesh.boundarycondition = [1;1;1;1]; % Set boundary condition for each boundary
mesh.curvedboundary = [];
mesh.curvedboundaryexpr = {};
mesh.periodicexpr = {};
mesh.porder = porder;
mesh.elemtype = elemtype;
mesh.f = facenumbering(mesh.p,mesh.t,mesh.elemtype,mesh.boundaryexpr,mesh.periodicexpr);
mesh.dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,mesh.porder);

x = (mesh.dgnodes(:,1,:));
y = (mesh.dgnodes(:,2,:));
r = sqrt(x.^2+y.^2);
a1 = 20;
a2 = 100;
a = 0.25;
rho = 1 + a1*sech(a2*(r.^2-a^2));
drhodx = -(2*a1*a2*x.*sinh(a2*(- a^2 + x.^2 + y.^2)))./cosh(a2*(- a^2 + x.^2 + y.^2)).^2;
drhody = -(2*a1*a2*y.*sinh(a2*(- a^2 + x.^2 + y.^2)))./cosh(a2*(- a^2 + x.^2 + y.^2)).^2;

