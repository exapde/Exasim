function mesh = mkmesh_isoq3d2(porder, ns, dR)

%dR = 0.01;
mesh2d = mkmesh_isoq2d2(porder, dR);
[p, t, xdg] = rotate_mesh(mesh2d.p, mesh2d.t, ns, mesh2d.dgnodes);

npe2d = size(mesh2d.dist,1);
npe3d = size(xdg,1);
nphi = npe3d/npe2d;
if abs(nphi-round(nphi)) > 1e-12
    error("mkmesh_isoq3d2: inconsistent DG-node counts when lifting mesh2d.dist to 3D.");
end
nphi = round(nphi);
dist = repmat(mesh2d.dist, [1, nphi, 1, ns]);
dist = reshape(dist, [npe3d, 1, size(xdg,3)]);

mesh = mkmesh(p',t',porder,{'true'},1,1);
mesh.p = mesh.p';
mesh.t = mesh.t';
mesh.dgnodes = xdg;
mesh.dist = dist;

%[~, xu] = isoq();

ymin = min(mesh.p(2,:));
zmin = min(mesh.p(3,:));
xmax = max(mesh.p(1,:));
% ymax = max(mesh.p(2,:));
% zmax = max(mesh.p(3,:));
tol = 1e-6;
% profiletol = 1e-4;
% ru = xu(:,2) + dR;
                     %@(p) bnd_box(p, [-1e-3 xmax 0 yc+1e-3 0 yc+1e-3], "inside"), ...

% @(p) -1e-3<p(1,:) && p(1,:)<xmax+1e-3 && -1e-3<p(2,:) && p(2,:)<yc+1e-3 && -1e-3<p(3,:) && p(3,:)<yc+1e-3 , ...
                     
% @(p) abs(sqrt(p(2,:).^2 + p(3,:).^2) - interp1(xu(:,1), ru, p(1,:), 'linear', 'extrap')) < profiletol, ... 

mesh.boundaryexpr = {@(p) abs(p(2,:)-ymin)<tol, ...
                     @(p) abs(p(3,:)-zmin)<tol, ...
                     @(p) abs(p(1,:)-xmax)<tol, ...
                     @(p) abs(p(2,:).^2 + p(3,:).^2 - dR^2)<tol, ...
                     @(p) -1e-3<p(1,:) && p(1,:)<xmax+1e-3 && sqrt(p(2,:).^2 + p(3,:).^2) < 0.06 + dR, ...                     
                     @(p) abs(p(1,:))< 20 + 1e-6};
mesh.boundarycondition = [1, 2, 1, 3, 2, 2]; % Split boundary 5 without changing the previous BC id
mesh.f = facenumbering(mesh.p,mesh.t,1,mesh.boundaryexpr,[]);
mesh.periodicboundary = [];
mesh.periodicexpr = {};

% figure(1); clf; meshplot(mesh);
% axis on; axis equal; axis tight;

%colors = ['b', 'r', 'g', 'y', 'm'];
colors = lines(12);
figure(2); clf; boundaryplot(mesh,1,colors(1,:));
hold on;
for i = 2:6
  boundaryplot(mesh,i,colors(i,:));
end
