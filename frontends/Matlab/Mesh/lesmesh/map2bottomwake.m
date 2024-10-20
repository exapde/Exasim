function p = map2bottomwake(pp,xf,yf)
p = 0*pp;
nr = 0*pp;
lb = pp(:,1) < 0;
% % clf;
xf =[xf, yf];
%nn = size(xf,1)-1;
t = 0:length(xf(:,1))-1;
spf = [spline(t,xf(:,1)),spline(t,xf(:,2))];
spfder = [fnder(spf(1)),fnder(spf(2))];
dxf = [fnval(spfder(1),t)', fnval(spfder(2),t)'];
dsf = sqrt(dxf(:,1).^2+dxf(:,2).^2);
dxf = [dxf(:,1)./dsf, dxf(:,2)./dsf];
nf = -[dxf(:,2), -dxf(:,1)];

%BOTTOM WAKE
tw = 0:6;
xwb = zeros(7,2);
xwb(1,:) = xf(1,:);
xwb(2,:) = [xf(1,1) + 0.5, xf(1,2)-0.01];
xwb(3,:) = [xf(1,1) + 1.5, xf(1,2)+0.0];
xwb(4,:) = [xf(1,1) + 3.0, xf(1,2)+0.04*1.5*1.1];
xwb(5,:) = [xf(1,1) + 5.0, xf(1,2)+0.04*3.5*1.25];
xwb(6,:) = [xf(1,1) + 7.5, xf(1,2)+0.04*6*1.5];
xwb(7,:) = [xf(1,1) + 10.0, xf(1,2)+0.04*8.5*1.7];
spwb = [spline(tw,xwb(:,1)),spline(tw,xwb(:,2))];
spwbder = [fnder(spwb(1)),fnder(spwb(2))];
% twp = 0:0.1:6;
%dxwb = [fnval(spwbder(1),twp)', fnval(spwbder(2),twp)'];
%dswb = sqrt(dxwb(:,1).^2+dxwb(:,2).^2);
%nb = [dxwb(:,2)./dswb, -dxwb(:,1)./dswb];
%nb = 0.5*(nb.*twp' + nf(1,:).*(2-twp'));
%nbs = sqrt(nb(:,1).^2+nb(:,2).^2);
%nb = 0.1*nb./nbs;

blth = 16;
p(lb,:) = [fnval(spwb(1),-6*pp(lb,1)),fnval(spwb(2),-6*pp(lb,1))];
nr(lb,:) = [fnval(spwbder(1),-6*pp(lb,1)), fnval(spwbder(2),-6*pp(lb,1))];
nrs = sqrt(nr(lb,1).^2+nr(lb,2).^2);
nr(lb,:) = [nr(lb,2)./nrs, -nr(lb,1)./nrs];
nr(lb,:) = (-nr(lb,:).*pp(lb,1) + nf(1,:).*(1+pp(lb,1)));
nrs = sqrt(nr(lb,1).^2+nr(lb,2).^2);
nr(lb,:) = nr(lb,:)./nrs;
p(lb,:) = p(lb,:) + (1+pp(lb,1)-pp(lb,1)*blth).*nr(lb,:).*pp(lb,2);

