function [xl, xu, Rn] = isoq()

Rn = 0.102;
Rs = Rn/16; 
xE = 0.06;
ThN = 27.82*pi/180;
Thf = pi/2;
n1 = 480;
n2 = 48;
n3 = 240;

th = 0:ThN/n1:ThN;
xc0 = Rn*(1 - cos(th));
yc0 = Rn*sin(th);

Thi = pi - ThN;
c1_x =  xc0(end) - Rs*cos(Thi);
c1_y =  yc0(end) - Rs*sin(Thi);

th = Thi:(Thf-Thi)/n2:Thf;
xc1_1 = Rs*cos(th);
yc1_1 = Rs*sin(th);

xc1_1 = xc1_1 + c1_x;
yc1_1 = yc1_1 + c1_y;

xl0 = [xc0(:) yc0(:)];
xl1 = [xc1_1(:) yc1_1(:)];

xT = max(xl1(:,1));
yT = max(xl1(:,2));

xl2 = [linspace(xT, xE, n3)' yT*ones(n3,1)];

figure(1); clf;
hold on;
plot(xl0(:,1), xl0(:,2), 'or');
plot(xl1(:,1), xl1(:,2), 'ob');
plot(xl2(:,1), xl2(:,2), 'og');
axis equal;

xl = [xl0(1:end-1,:); xl1(1:end-1,:); xl2(1:end,:)];

t = linspace(pi, pi/2, 2000)';
xc = 3.8*xE + 5*xE*cos(t);
yc = 3.5*xE*sin(t);
xu = [xc yc];

ind = xu(:,1) <= xE;
xu = xu(ind,:);
xu(end,1) = xE;

figure(2); clf;
hold on;
plot(xl(:,1), xl(:,2), 'or');
plot(xu(:,1), xu(:,2), 'ob');
axis equal;



% 
% plot(xl2(:,1), xl2(:,2), 'og');
% plot(xl3(:,1), xl3(:,2), 'ok');
% 
% plot(c1_x, c1_y, '-ok');
% 
% %xo = [xl0(1:end-1,:); xl1(1:end-1,:); xl2(1:end-1,:); xl3(1:end,:)];
% 
% n1 = 5000;
% n2 = 5000;
% n3 = 2000;
% n4 = 1000;
% 
% ThN = 23.0353*pi/180;
% th = 0:ThN/n1:ThN;
% xc0 = Rn*(1 - cos(th));
% yc0 = Rn*sin(th);
% xr0 = [xc0(:) yc0(:)];
% 
% Thi = pi-ThN;
% Thf = pi/2-32.5*pi/180;
% th = Thi:(Thf-Thi)/n2:Thf;
% xc1 = Rs*cos(th) + c1_x;
% yc1 = Rs*sin(th) + c1_y;
% xr1 = [xc1(:) yc1(:)];
% 
% xi = linspace(0,1,n3)';
% XA = xr1(end,:);
% XB = [3.3020 0.9116];
% x = XA(1,1) + (XB(1,1)-XA(1,1))*xi;
% y = XA(1,2) + (XB(1,2)-XA(1,2))*xi;
% xr2 = [x y];
% 
% xi = linspace(0,1,n4)';
% XA = XB;
% XB = [3.3020 0];
% x = XA(1,1) + (XB(1,1)-XA(1,1))*xi;
% y = XA(1,2) + (XB(1,2)-XA(1,2))*xi;
% xr3 = [x y];
% 
% % figure(2); clf;
% % hold on;
% % plot(xr0(:,1), xr0(:,2), '-b');
% % plot(xr1(:,1), xr1(:,2), '-b');
% % plot(xr2(:,1), xr2(:,2), '-b');
% % plot(xr3(:,1), xr3(:,2), '-b');
% % plot(c1_x, c1_y, '-ok');
% 
% xr = [xr0(1:end-1,:); xr1(1:end-1,:); xr2(1:end,:)];
% 
% xs = shifted_curve_normal_direction(xr, D);
% xs(1,:) = [D 0];
% XC = XB;
% XC(1) = XB(1) - D;
% s = (xs(end-3,2)-xs(end-4,2))/(xs(end-3,1)-xs(end-4,1));
% c =  xs(end-4,2) - s*xs(end-4,1);
% XC(2) = c + s * XC(1);
% 
% % [xr(1:10,:) xs(1:10,:)]
% % pause
% 
% ind = xs(:,1) < XC(1);
% xs = xs(ind,:);
% 
% xi = linspace(0,1,n4)';
% XA = XC;
% XB(1) = XB(1) - D;
% x = XA(1,1) + (XB(1,1)-XA(1,1))*xi;
% y = XA(1,2) + (XB(1,2)-XA(1,2))*xi;
% xs3 = [x y];
% xs = [xs; xs3];
% 
% xf = [xr0(1:end-1,:); xr1(1:end-1,:); xr2(1:end-1,:); xr3(1:end,:)];
% 
% figure(2); clf;
% hold on;
% plot(xf(:,1), xf(:,2), '-b');
% plot(xs(:,1), xs(:,2), '-r');
% %plot(XC(1), XC(2), 'ok');
% axis equal; axis tight;
% 
