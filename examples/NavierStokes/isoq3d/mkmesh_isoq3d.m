function mesh = mkmesh_isoq3d(porder, ns)

[mesh2d, xc, xd] = mkmesh_isoq2d(porder);
[p, t, xdg] = rotate_mesh(mesh2d.p, mesh2d.t, ns, mesh2d.dgnodes);

mesh = mkmesh(p',t',porder,{'true'},1,1);
mesh.p = mesh.p';
mesh.t = mesh.t';
mesh.dgnodes = xdg;

xc = [xc; 0*xc(1,:)];
xd = [xd; 0*xd(1,:)];
[xl, ~, Rn, xc0, yc0] = isoq();
mesh2 = mkmesh_sphericalfrustum(xc, porder, ns, Rn, xd);

figure(1); clf; meshplot(mesh2d);
axis on; axis equal; %axis tight;
hold on;
x = mesh2.dgnodes(:,1,:);
y = mesh2.dgnodes(:,2,:);
z = mesh2.dgnodes(:,3,:);
ind = abs(z(:)) <= 1e-8;
plot(x(ind), y(ind), 'ok');
x = mesh2d.dgnodes(:,1,:);
y = mesh2d.dgnodes(:,2,:);
plot(x(:), y(:), 'or');

% figure(1); clf; meshplot(mesh);
% axis on; axis equal; %axis tight;
% hold on;
% meshplot(mesh2);
% pause

[p,t] = connectmesh(mesh2.p',mesh2.t',mesh.p',mesh.t');
mesh = mkmesh(p,t,porder,{'true'},1,1);
mesh.p = mesh.p';
mesh.t = mesh.t';
mesh.dgnodes = cat(3, mesh2.dgnodes, xdg);

% figure(2); clf; meshplot(mesh,1);
% axis on; axis equal; axis tight;

% figure(3); clf; meshplot(mesh2d,1);
% axis on; axis equal; axis tight; hold on;
% plot(xc(1,:),xc(2,:),'o');
% plot(xd(1,:),xd(2,:),'*');
% x = mesh2.dgnodes(:,1,:);
% y = sqrt(mesh2.dgnodes(:,2,:).^2 + mesh2.dgnodes(:,3,:).^2);
% plot(x(:),y(:),'s');
% plot(xc0(:),yc0(:),'-r', 'LineWidth', 2);


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
mesh.boundarycondition = [1, 2, 1, 3, 2]; % Set boundary condition for each boundary
mesh.f = facenumbering(mesh.p,mesh.t,1,mesh.boundaryexpr,[]);
mesh.periodicboundary = [];
mesh.periodicexpr = {};

% figure(1); clf; meshplot(mesh,1);
% axis on; axis equal; axis tight;

%colors = ['b', 'r', 'g', 'y', 'm'];
% colors = lines(12);
% figure(2); clf; boundaryplot(mesh,1,colors(1,:));
% hold on;
% for i = 2:5
%   boundaryplot(mesh,i,colors(i,:));
% end
