function p = map2foil(pp,xf,yf)

p = 0*pp;
nr = 0*pp;
lf = 1:size(pp,1);
xf =[xf, yf];
nn = size(xf,1)-1;
t = 0:length(xf(:,1))-1;
spf = [spline(t,xf(:,1)),spline(t,xf(:,2))];
spfder = [fnder(spf(1)),fnder(spf(2))];
p(lf,1:2) = [fnval(spf(1),nn*pp(lf,1)/2),fnval(spf(2),nn*pp(lf,1)/2)];
nr(lf,1:2) = [fnval(spfder(1),nn*pp(lf,1)/2), fnval(spfder(2),nn*pp(lf,1)/2)];
nrs = sqrt(nr(lf,1).^2+nr(lf,2).^2);
nr(lf,1:2) = -[nr(lf,2)./nrs, -nr(lf,1)./nrs];
p(lf,1:2) = p(lf,1:2) + nr(lf,1:2).*pp(lf,2);

