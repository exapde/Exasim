function mesh = mkmesh_isoq3d(porder, ns)


[mesh2d, xc, xd] = mkmesh_isoq2d(porder);

xc = [xc; 0*xc(1,:)];
xd = [xd; 0*xd(1,:)];
[xl, ~, Rn] = isoq();
mesh1 = mkmesh_sphericalfrustum(xc, porder, ns, Rn, 0, xd);

[mesh2.p, mesh2.t, mesh2.dgnodes] = rotate_mesh(mesh2d.p, mesh2d.t, ns, mesh2d.dgnodes);
% mesh2 = mkmesh(mesh2.p',mesh2.t',porder,{'true'},1,1);
% mesh2.p = mesh2.p';
% mesh2.t = mesh2.t';

% figure(1); clf; meshplot(mesh2,1);
% axis on; axis equal; %axis tight;
% hold on;
% meshplot(mesh1,2);

[p,t] = connectmesh(mesh1.p',mesh1.t',mesh2.p',mesh2.t');
mesh = mkmesh(p,t,porder,{'true'},1,1);
mesh.p = mesh.p';
mesh.t = mesh.t';
mesh.dgnodes = cat(3, mesh1.dgnodes, mesh2.dgnodes);

%xmin = min(mesh.p(1,:));
ymin = min(mesh.p(2,:));
zmin = min(mesh.p(3,:));

xmax = max(mesh.p(1,:));
% ymax = max(mesh.p(2,:));
% zmax = max(mesh.p(3,:));
yc =  max(xl(:,2));

                     %@(p) bnd_box(p, [-1e-3 xmax 0 yc+1e-3 0 yc+1e-3], "inside"), ...
mesh.boundaryexpr = {@(p) abs(p(2,:)-ymin)<1e-6, ...
                     @(p) abs(p(3,:)-zmin)<1e-6, ...
                     @(p) abs(p(1,:)-xmax)<1e-6, ...
                     @(p) -1e-3<p(1,:) && p(1,:)<xmax+1e-3 && -1e-3<p(2,:) && p(2,:)<yc+1e-3 && -1e-3<p(3,:) && p(3,:)<yc+1e-3 , ...
                     @(p) abs(p(1,:))< 20 + 1e-6};
% axis symmetric, outflow, inflow, wall
mesh.boundarycondition = [1, 2, 1, 3, 2]; % Set boundary condition for each boundary
mesh.f = facenumbering(mesh.p,mesh.t,1,mesh.boundaryexpr,[]);
mesh.periodicboundary = [];
mesh.periodicexpr = {};

figure(1); clf; meshplot(mesh,1);
axis on; axis equal; axis tight;

%colors = ['b', 'r', 'g', 'y', 'm'];
colors = lines(12);
figure(2); clf; boundaryplot(mesh,1,colors(1,:));
hold on;
for i = 2:5
  boundaryplot(mesh,i,colors(i,:));
end
