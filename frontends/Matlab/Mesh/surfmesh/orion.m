load orioncurves.mat
load xm.mat;
xu = xm;

xmax = 3.2;
x1 = -0.2;
x2 = 0.85;
x3 = 0.4812+0.0187;
x4 = 0.8477+0.0256;

ind = (xm(:,1) <= xmax);
xm = xm(ind,:);

ind = (xl(:,1) <= xmax);
xl = xl(ind,:);

ind = (xu(:,1) <= xmax);
xu = xu(ind,:);

%figure(1);clf;plot(xl(:,1),xl(:,2),xu(:,1),xu(:,2));

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

figure(1);clf;plot(xu1(:,1),xu1(:,2),xu2(:,1),xu2(:,2),xu3(:,1),xu3(:,2));
figure(2);clf;plot(xl1(:,1),xl1(:,2),xl2(:,1),xl2(:,2),xl3(:,1),xl3(:,2));

[p1, t, dgnodes1, mesh1, il1, iu1] = surfmesh2d(xl1, xu1, 48, 60, 3, [2 2], [9, 6]);
[p2, t, dgnodes2, mesh2, il2, iu2] = surfmesh2d(xl2, xu2, 24, 60, 3, [1 1], [9, 6], 0);
[p3, t, dgnodes3, mesh3, il3, iu3] = surfmesh2d(xl3, xu3, 24, 60, 3, [3 1], [9, 6], 0.6);

pu1 = p1(:,iu1)';
pu2 = p2(:,iu2)';
pu3 = p3(:,iu3)';

figure(3); clf; meshplot(mesh1,0);
hold on;
meshplot(mesh2,0);
meshplot(mesh3,0);

% plot(pu1(:,1), pu1(:,2), 'o', 'LineWidth', 2);
% plot(pu2(:,1), pu2(:,2), 'o', 'LineWidth', 2);
% plot(pu3(:,1), pu3(:,2), 'o', 'LineWidth', 2);

% plot(xl1(:,1), xl1(:,2), 'o', 'LineWidth', 2);
% plot(xu1(:,1), xu1(:,2), 's', 'LineWidth', 2);
% plot(xl2(:,1), xl2(:,2), 'o', 'LineWidth', 2);
% plot(xu2(:,1), xu2(:,2), 's', 'LineWidth', 2);
% plot(xl3(:,1), xl3(:,2), 'o', 'LineWidth', 2);
% plot(xu3(:,1), xu3(:,2), 's', 'LineWidth', 2);

ymin = 1e-4;
x = mesh3.dgnodes(:,1,:);
xmax = max(x(:));
t = linspace(pi, pi/2, 1000);
xc = 5.8 + 7*cos(t);
yc = ymin + 6.8*sin(t);
ind = xc <=xmax;
xc = xc(ind);
yc = yc(ind);

% plot(xm(:,1), xm(:,2), '-', 'LineWidth', 3);
% plot(xc, yc, '-');

axis equal; axis tight;
% 
% pu = unique([pu1; pu2; pu3], 'rows');
% [p, t, dgnodes4, mesh4] = surfmesh2d(xm, [xc(:) yc(:)], pu, 10, 3, [0 0], [5, 1]);
% %figure(2); clf; 
% meshplot(mesh4,0);
% %hold on;
% %plot(pu(:,1), pu(:,2), 'o', 'LineWidth', 2);
% x = squeeze(dgnodes3(:,1,:));
% y = squeeze(dgnodes3(:,2,:));
% plot(x(:),y(:),'o');
% x = squeeze(dgnodes4(:,1,:));
% y = squeeze(dgnodes4(:,2,:));
% plot(x(:),y(:),'o');



