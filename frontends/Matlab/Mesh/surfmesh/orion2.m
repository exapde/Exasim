load orioncurves.mat
load xm.mat;

ymin = 1e-4;
xmax = 3.2;
x1 = -0.25;
x2 = 0.5;
x3 = 0.4812+0.0187;
x4 = 0.8477+0.0256;

t = linspace(pi, pi/2, 1000)';
xc = 5.8 + 7*cos(t);
yc = ymin + 6.8*sin(t);
ind = xc <=xmax;
xu = [xc(ind) yc(ind)];

ind = (xm(:,1) <= xmax);
xm = xm(ind,:);

ind = (xl(:,1) <= xmax);
xl = xl(ind,:);

ind = (xu(:,1) <= xmax);
xu = xu(ind,:);

figure(1);clf;plot(xl(:,1),xl(:,2),xu(:,1),xu(:,2),xm(:,1),xm(:,2));

ind = (xu(:,1) <= x1);
xu1 = xu(ind,:);

ind = (xu(:,1) > x1 & xu(:,1) <= x2);
xu2 = [xu1(end,:); xu(ind,:)];

ind = (xu(:,1) > x2);
xu3 = [xu2(end,:); xu(ind,:)];

ind = (xl(:,1) <= x3);
xl1 = xl(ind,:);

ind = (xl(:,1) > x3 & xl(:,1) <= x4);
xl2 = [xl1(end,:); xl(ind,:)];

ind = (xl(:,1) > x4);
xl3 = [xl2(end,:); xl(ind,:)];

xm1 = [0 0];
j = 1;
for i = 1:size(xm,1)
  [~, ~, position] = find_perpendicular_point(xl1(end,1), xl1(end,2), xu1(end,1), xu1(end,2), xm(i,1), xm(i,2));
  if position<=0
    xm1(j,:) = xm(i,:);
    j = j + 1;
  end
end

xm2 = [0 0];
j = 1;
for i = 1:size(xm,1)
  [~, ~, position1] = find_perpendicular_point(xl1(end,1), xl1(end,2), xu1(end,1), xu1(end,2), xm(i,1), xm(i,2));  
  [~, ~, position2] = find_perpendicular_point(xl2(end,1), xl2(end,2), xu2(end,1), xu2(end,2), xm(i,1), xm(i,2));
  if (position1>=0) && (position2<=0)
    xm2(j,:) = xm(i,:);
    j = j + 1;
  end
end

xm3 = [0 0];
j = 1;
for i = 1:size(xm,1)
  [~, ~, position] = find_perpendicular_point(xl2(end,1), xl2(end,2), xu2(end,1), xu2(end,2), xm(i,1), xm(i,2));
  if position>=0
    xm3(j,:) = xm(i,:);
    j = j + 1;
  end
end

mesh1 = surfmesh2d_shockcurve(xl1, xu1, xm1, 48, 80, 10, 3, [2.5 2], [8 6.4], [4 0], 0);
mesh2 = surfmesh2d_shockcurve(xl2, xu2, xm2, 20, 80, 10, 3, [.1 .1], [8 6.4], [4 0], 0);
mesh3 = surfmesh2d_shockcurve(xl3, xu3, xm3, 24, 80, 10, 3, [3 1], [8 6.4], [4 0], 0.75);

% figure(1);clf;plot(xl1(:,1),xl1(:,2),xu1(:,1),xu1(:,2),xm1(:,1),xm1(:,2));
% figure(3);clf;plot(xl2(:,1),xl2(:,2),xu2(:,1),xu2(:,2),xm2(:,1),xm2(:,2));

figure(1); clf; 
meshplot(mesh1);
hold on;
meshplot(mesh2);
meshplot(mesh3);
% x=mesh1.dgnodes(:,1,:); y=mesh1.dgnodes(:,2,:); plot(x(:),y(:),'o','linewidth',1);
% x=mesh2.dgnodes(:,1,:); y=mesh2.dgnodes(:,2,:); plot(x(:),y(:),'s','linewidth',1);
% x=mesh3.dgnodes(:,1,:); y=mesh3.dgnodes(:,2,:); plot(x(:),y(:),'d','linewidth',1);
plot(xm(:,1),xm(:,2),'-r','linewidth',2);
axis equal; axis tight;

% mesh.xpe = mesh.plocal;
% figure(1); clf; meshplot(mesh);

mesh = mesh1;
[mesh.p, mesh.t] = connectmesh(mesh1.p', mesh1.t', mesh2.p', mesh2.t', 1e-6);
[mesh.p, mesh.t] = connectmesh(mesh.p, mesh.t, mesh3.p', mesh3.t', 1e-6);
mesh.dgnodes = cat(3, mesh1.dgnodes, mesh2.dgnodes);
mesh.dgnodes = cat(3, mesh.dgnodes, mesh3.dgnodes);

mesh.p = mesh.p';
mesh.t = mesh.t';
deltay = min(mesh.p(2,:));
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(2,:)-deltay)<1e-4, @(p) p(1,:)> 3.125+1e-6, @(p) ((p(1,:) < 1e-3) | (p(2,:) > 2.6)), @(p) abs(p(1,:))< 10 + 1e-6};
% axis symmetric, outflow, inflow, wall
mesh.boundarycondition = [6, 2, 1, 3]; % Set boundary condition for each boundary
mesh.f = facenumbering(mesh.p,mesh.t,1,mesh.boundaryexpr,[]);

% figure(1);clf;
% hold on;
% boundaryplot(mesh,1);
% boundaryplot(mesh,2);
% boundaryplot(mesh,3);
% boundaryplot(mesh,4);
% for i = 1:4
%   boundaryplot(mesh,i);
% end


