function [xl, xs, xf] = orion(D)

if nargin<1
  D = 0.2;
end

Rn = 6.0350;
Rs = 0.2515;

ThN = 23.0353*pi/180;
th = 0:ThN/480:ThN;
xc0 = Rn*(1 - cos(th));
yc0 = Rn*sin(th);

Thi = pi - ThN;
c1_x =  xc0(end) - Rs*cos(Thi);
c1_y =  yc0(end) - Rs*sin(Thi);

XA = -0.1;   % Outer boudnary location
YA =  4.0;
Thf = atan2(YA - c1_y, XA - c1_x);

th = Thi:(Thf-Thi)/40:Thf;
xc1_1 = Rs*cos(th);
yc1_1 = Rs*sin(th);

xc1_1 = xc1_1 + c1_x;
yc1_1 = yc1_1 + c1_y;


Thi = Thf;
ThA = 32.5*pi/180;
Thf = pi/2 - ThA;

th = Thi:(Thf-Thi)/40:Thf;
xc1_2 = Rs*cos(th);
yc1_2 = Rs*sin(th);

xc1_2 = xc1_2 + c1_x;
yc1_2 = yc1_2 + c1_y;

xl0 = [xc0(:) yc0(:)];
xl1 = [xc1_1(:) yc1_1(:); xc1_2(:) yc1_2(:)];

xi = linspace(0,1,200)';
XA = xl1(end,:);
%XB = [3.2785 0.9266];
XB = [3.3020 0.9116];
x = XA(1,1) + (XB(1,1)-XA(1,1))*xi;
y = XA(1,2) + (XB(1,2)-XA(1,2))*xi;
xl2 = [x y];

XA = XB;
XB = [3.3020 0];

xi = linspace(0,1,50)';
x = XA(1,1) + (XB(1,1)-XA(1,1))*xi;
y = XA(1,2) + (XB(1,2)-XA(1,2))*xi;
xl3 = [x y];

figure(1); clf;
hold on;
plot(xl0(:,1), xl0(:,2), 'or');
plot(xl1(:,1), xl1(:,2), 'ob');
plot(xl2(:,1), xl2(:,2), 'og');
plot(xl3(:,1), xl3(:,2), 'ok');

plot(c1_x, c1_y, '-ok');

xl = [xl0(1:end-1,:); xl1(1:end-1,:); xl2(1:end,:)];
%xo = [xl0(1:end-1,:); xl1(1:end-1,:); xl2(1:end-1,:); xl3(1:end,:)];

n1 = 5000;
n2 = 5000;
n3 = 2000;
n4 = 1000;

ThN = 23.0353*pi/180;
th = 0:ThN/n1:ThN;
xc0 = Rn*(1 - cos(th));
yc0 = Rn*sin(th);
xr0 = [xc0(:) yc0(:)];

Thi = pi-ThN;
Thf = pi/2-32.5*pi/180;
th = Thi:(Thf-Thi)/n2:Thf;
xc1 = Rs*cos(th) + c1_x;
yc1 = Rs*sin(th) + c1_y;
xr1 = [xc1(:) yc1(:)];

xi = linspace(0,1,n3)';
XA = xr1(end,:);
XB = [3.3020 0.9116];
x = XA(1,1) + (XB(1,1)-XA(1,1))*xi;
y = XA(1,2) + (XB(1,2)-XA(1,2))*xi;
xr2 = [x y];

xi = linspace(0,1,n4)';
XA = XB;
XB = [3.3020 0];
x = XA(1,1) + (XB(1,1)-XA(1,1))*xi;
y = XA(1,2) + (XB(1,2)-XA(1,2))*xi;
xr3 = [x y];

% figure(2); clf;
% hold on;
% plot(xr0(:,1), xr0(:,2), '-b');
% plot(xr1(:,1), xr1(:,2), '-b');
% plot(xr2(:,1), xr2(:,2), '-b');
% plot(xr3(:,1), xr3(:,2), '-b');
% plot(c1_x, c1_y, '-ok');

xr = [xr0(1:end-1,:); xr1(1:end-1,:); xr2(1:end,:)];

xs = shifted_curve_normal_direction(xr, D);
xs(1,:) = [D 0];
XC = XB;
XC(1) = XB(1) - D;
s = (xs(end-3,2)-xs(end-4,2))/(xs(end-3,1)-xs(end-4,1));
c =  xs(end-4,2) - s*xs(end-4,1);
XC(2) = c + s * XC(1);

% [xr(1:10,:) xs(1:10,:)]
% pause

ind = xs(:,1) < XC(1);
xs = xs(ind,:);

xi = linspace(0,1,n4)';
XA = XC;
XB(1) = XB(1) - D;
x = XA(1,1) + (XB(1,1)-XA(1,1))*xi;
y = XA(1,2) + (XB(1,2)-XA(1,2))*xi;
xs3 = [x y];
xs = [xs; xs3];

xf = [xr0(1:end-1,:); xr1(1:end-1,:); xr2(1:end-1,:); xr3(1:end,:)];

figure(2); clf;
hold on;
plot(xf(:,1), xf(:,2), '-b');
plot(xs(:,1), xs(:,2), '-r');
%plot(XC(1), XC(2), 'ok');
axis equal; axis tight;

