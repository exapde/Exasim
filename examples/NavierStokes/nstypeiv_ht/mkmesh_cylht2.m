function mesh = mkmesh_cylht(porder)

mesh = mkmesh_square(4,60,porder,1,1,1,1,1);
% mesh.p(1,:) = loginc(mesh.p(1,:), 3);
% mesh.dgnodes(:,1,:) = loginc(mesh.dgnodes(:,1,:), 3);
mesh = mkmesh_halfcircle(mesh, 0.5, 1.0, 1.0, pi/2, 3*pi/2);
mesh.porder = porder;
mesh.boundaryexpr = {@(p) sqrt(p(1,:).^2+p(2,:).^2)>1-1e-3, @(p) p(1,:)>-1e-3, @(p) abs(p(1,:))<20};
mesh.periodicexpr = {};
