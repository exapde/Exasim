function p = map(pp)
p = 0*pp;
nr = 0*pp;
lf = (pp(:,1) >=0 & pp(:,1) <=2);
lb = pp(:,1) < 0;
lt = pp(:,1) > 2;
load('foilgeom.mat')
% % clf;
xf =[xf, yf];
nn = size(xf,1)-1;
clear yf;
t = 0:length(xf(:,1))-1;
spf = [spline(t,xf(:,1)),spline(t,xf(:,2))];
spfder = [fnder(spf(1)),fnder(spf(2))];
dxf = [fnval(spfder(1),t)', fnval(spfder(2),t)'];
dsf = sqrt(dxf(:,1).^2+dxf(:,2).^2);
dxf = [dxf(:,1)./dsf, dxf(:,2)./dsf];
nf = -[dxf(:,2), -dxf(:,1)];
p(lf,:) = [fnval(spf(1),nn*pp(lf,1)/2),fnval(spf(2),nn*pp(lf,1)/2)];
nr(lf,:) = [fnval(spfder(1),nn*pp(lf,1)/2), fnval(spfder(2),nn*pp(lf,1)/2)];
nrs = sqrt(nr(lf,1).^2+nr(lf,2).^2);
nr(lf,:) = -[nr(lf,2)./nrs, -nr(lf,1)./nrs];
p(lf,:) = p(lf,:) + nr(lf,:).*pp(lf,2);
% plot(xf(:,1),xf(:,2),'-o',xf(:,1)+0.1*nf(:,1),xf(:,2)+0.1*nf(:,2),'-o');
% hold on;
% axis equal;
%BOTTOM WAKE
tw = 0:6;
xwb = zeros(7,2);
xwb(1,:) = xf(1,:);
xwb(2,:) = [xf(1,1) + 0.5, xf(1,2)-0.01];
xwb(3,:) = [xf(1,1) + 1.5, xf(1,2)+0.0];
xwb(4,:) = [xf(1,1) + 3.0, xf(1,2)+0.04*1.5];
xwb(5,:) = [xf(1,1) + 5.0, xf(1,2)+0.04*3.5];
xwb(6,:) = [xf(1,1) + 7.5, xf(1,2)+0.04*6];
xwb(7,:) = [xf(1,1) + 10.0, xf(1,2)+0.04*8.5];
spwb = [spline(tw,xwb(:,1)),spline(tw,xwb(:,2))];
spwbder = [fnder(spwb(1)),fnder(spwb(2))];
twp = 0:0.1:6;
xwbp = [fnval(spwb(1),twp)', fnval(spwb(2),twp)'];
dxwb = [fnval(spwbder(1),twp)', fnval(spwbder(2),twp)'];
dswb = sqrt(dxwb(:,1).^2+dxwb(:,2).^2);
nb = [dxwb(:,2)./dswb, -dxwb(:,1)./dswb];
nb = 0.5*(nb.*twp' + nf(1,:).*(2-twp'));
nbs = sqrt(nb(:,1).^2+nb(:,2).^2);
nb = 0.1*nb./nbs;
% plot(xwbp(:,1),xwbp(:,2),'-+',xwbp(:,1)+nb(:,1),xwbp(:,2)+nb(:,2),'-+')
%
% for i = 1:length(twp)
% vectarrow([xwbp(i,1),xwbp(i,2)],[xwbp(i,1)+nb(i,1),xwbp(i,2)+nb(i,2)]);
% hold on;
% end
blth = 6;
p(lb,:) = [fnval(spwb(1),-6*pp(lb,1)),fnval(spwb(2),-6*pp(lb,1))];
nr(lb,:) = [fnval(spwbder(1),-6*pp(lb,1)), fnval(spwbder(2),-6*pp(lb,1))];
nrs = sqrt(nr(lb,1).^2+nr(lb,2).^2);
nr(lb,:) = [nr(lb,2)./nrs, -nr(lb,1)./nrs];
nr(lb,:) = (-nr(lb,:).*pp(lb,1) + nf(1,:).*(1+pp(lb,1)));
nrs = sqrt(nr(lb,1).^2+nr(lb,2).^2);
nr(lb,:) = nr(lb,:)./nrs;
p(lb,:) = p(lb,:) + (1+pp(lb,1)-pp(lb,1)*blth).*nr(lb,:).*pp(lb,2);
%TOP WAKE
xwt = zeros(7,2);
xwt(1,:) = xf(end,:);
xwt(2,:) = [xf(end,1) + 0.5, xf(end,2)-0.0];
xwt(3,:) = [xf(end,1) + 1.5, xf(end,2)+0.03];
xwt(4,:) = [xf(end,1) + 3.0, xf(end,2)+0.03+0.06*1.5];
xwt(5,:) = [xf(end,1) + 5.0, xf(end,2)+0.03+0.06*3.5];
xwt(6,:) = [xf(end,1) + 7.5, xf(end,2)+0.03+0.06*6];
xwt(7,:) = [xf(end,1) + 10.0,xf(end,2)+0.03+0.06*8.5];
spwt = [spline(tw,xwt(:,1)),spline(tw,xwt(:,2))];
spwtder = [fnder(spwt(1)),fnder(spwt(2))];
xwtp = [fnval(spwt(1),twp)', fnval(spwt(2),twp)'];
dxwt = [fnval(spwtder(1),twp)', fnval(spwtder(2),twp)'];
dswt = sqrt(dxwt(:,1).^2+dxwt(:,2).^2);
nt = -[dxwt(:,2)./dswt, -dxwt(:,1)./dswt];
nt = 0.5*(nt.*twp' + nf(end,:).*(2-twp'));
nts = sqrt(nt(:,1).^2+nt(:,2).^2);
nt = 0.1*nt./nts;
p(lt,:) = [fnval(spwt(1),6*(pp(lt,1)-2)),fnval(spwt(2),6*(pp(lt,1)-2))];
nr(lt,:) = [fnval(spwtder(1),6*(pp(lt,1)-2)), fnval(spwtder(2),6*(pp(lt,1)-2))];
nrs = sqrt(nr(lt,1).^2+nr(lt,2).^2);
nr(lt,:) = -[nr(lt,2)./nrs, -nr(lt,1)./nrs];
nr(lt,:) = (nr(lt,:).*(pp(lt,1)-2) + nf(end,:).*(3-pp(lt,1)));
nrs = sqrt(nr(lt,1).^2+nr(lt,2).^2);
nr(lt,:) = nr(lt,:)./nrs;
p(lt,:) = p(lt,:) + ((3-pp(lt,1))+(pp(lt,1)-2)*blth).*nr(lt,:).*pp(lt,2);
% plot(xwtp(:,1),xwtp(:,2),'-+',xwtp(:,1)+nt(:,1),xwtp(:,2)+nt(:,2),'-+')
% for i = 1:length(twp)
% vectarrow([xwtp(i,1),xwtp(i,2)],[xwtp(i,1)+nt(i,1),xwtp(i,2)+nt(i,2)]);
% hold on;
% end
% nw = 20;
% ttwp = 0:20;
% d0 = 0.025;
% rat = 1.05;
% ttwp = d0*(rat.^ttwp -1)/(rat-1);
% ttwp = 2*ttwp/ttwp(end);
%
% xwttp = [fnval(spwt(1),ttwp)', fnval(spwt(2),ttwp)'];
% plot(xwttp(:,1),xwttp(:,2),'*');
% xwtbp = [fnval(spwb(1),ttwp)', fnval(spwb(2),ttwp)'];
% plot(xwtbp(:,1),xwtbp(:,2),'*');
%
%
% axis equal;
% grid on;
% xwt = xf(1,:) - xw(1,:);
% xwb = xf(end,:) - xw(1,:);
%
% spwx = spline(tw,xw(:,1));
% spwy = spline(tw,xw(:,2));
% spwxder = fnder(spwx);
% spwyder = fnder(spwy);
%
% xwp = zeros(length(ttwp),2);
% xwp(:,1) = fnval(spwx,ttwp);
% xwp(:,2) = fnval(spwy,ttwp);
% dwx = fnval(spwxder,ttwp);
% dwy = fnval(spwyder,ttwp);
% dws = sqrt(dwx.^2 + dwy.^2);
% dwx = dwx./dws;
% dwy = dwy./dws;
% xwt = xwp + xwt;
% xwb = xwp + xwb;
% xwot = xwt + 0.1*[-dwy' + dwx'];
% xwob = xwt - 0.1*[-dwy' + dwx'];
% plot(xwot(:,1),xwot(:,2),'-o',xwob(:,1),xwob(:,2),'-o');
%
%
%
%
% ttp = distribute(80,spx,spy,n);
% xp = fnval(spx,ttp);
% yp = fnval(spy,ttp);
%
% plot(xp,yp,'*');
% hold on;
% axis equal;
%
% dx = fnval(spxder,t);
% dy = fnval(spyder,t);
% ds = sqrt(dx.^2 + dy.^2);
% dx = dx./ds;
% dy = dy./ds;
% xfo = xf(:,1)-0.1*dy';
% yfo = xf(:,2)+0.1*dx';
%
% spxo = spline(t,xfo);
% spyo = spline(t,yfo);
% xpo = fnval(spxo,t);
% ypo = fnval(spyo,t);
% plot(xpo,ypo,'o');
%
% for i = 1:length(t)
% vectarrow([xf(i),yf(i)],[xfo(i),yfo(i)]);
% hold on;
% end
%
% axis equal;
