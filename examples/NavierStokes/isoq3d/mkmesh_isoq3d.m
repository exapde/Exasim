function mesh = mkmesh_isoq3d(porder, ns)


[mesh2d, xc, xd] = mkmesh_isoq2d(porder);
[p, t, xdg] = rotate_mesh(mesh2d.p, mesh2d.t, ns, mesh2d.dgnodes);

mesh = mkmesh(p',t',porder,{'true'},1,1);
mesh.p = mesh.p';
mesh.t = mesh.t';
mesh.dgnodes = xdg;

xc = [xc; 0*xc(1,:)];
xd = [xd; 0*xd(1,:)];
[~, ~, Rn] = isoq();
mesh2 = mkmesh_sphericalfrustum(xc, porder, ns, Rn, 0, xd);

figure(1); clf; meshplot(mesh,1);
axis on; axis equal; %axis tight;
hold on;
meshplot(mesh2,2);

[p,t] = connectmesh(mesh2.p',mesh2.t',mesh.p',mesh.t');
mesh = mkmesh(p,t,porder,{'true'},1,1);
mesh.p = mesh.p';
mesh.t = mesh.t';
mesh.dgnodes = cat(3, mesh2.dgnodes, xdg);

figure(2); clf; meshplot(mesh,1);
axis on; axis equal; axis tight;
