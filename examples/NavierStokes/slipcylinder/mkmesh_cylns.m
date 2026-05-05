function mesh = mkmesh_cylns(porder)

mesh = mkmesh_square(121,60,porder,1,1,1,1,1);
mesh.p(1,:) = logdec(mesh.p(1,:), 8);
mesh.dgnodes(:,1,:) = logdec(mesh.dgnodes(:,1,:), 8);
mesh = mkmesh_halfcircle(mesh, 1, 2.2, 3.5, pi/2, 3*pi/2);
mesh.porder = porder;
mesh.boundaryexpr = {@(p) sqrt(p(1,:).^2+p(2,:).^2)<1+1e-6, @(p) p(1,:)>-1e-7, @(p) abs(p(1,:))<20};
mesh.periodicexpr = {};

figure(1); clf; meshplot(mesh);