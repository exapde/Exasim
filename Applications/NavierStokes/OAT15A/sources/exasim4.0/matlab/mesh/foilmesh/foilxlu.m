function [xfl,yfl,xfu,yfu] = foilxlu(xf,yf,x1,y1,x2,y2)

% x1 = 0.7; y1 = 0.04;
% x2 = 1.0; y2 = 0.00;

c = [x1 1; x2 1]\[y1; y2];

xf = xf(:);
yf = yf(:);

ind = xf<x1;
xf1 = xf(ind);
yf1 = yf(ind);
[xf1,ind] = sort(xf1);
yf1 = yf1(ind); 

ind = yf1<=0;
xf1l = xf1(ind);
yf1l = yf1(ind);

ind = yf1>=0;
xf1u = xf1(ind);
yf1u = yf1(ind);

ind = xf>=x1;
xf2 = xf(ind);
yf2 = yf(ind);
[xf2,ind] = sort(xf2);
yf2 = yf2(ind); 

ind = ((c(1)*xf2 + c(2) - yf2)>0);
xf2l = xf2(ind);
yf2l = yf2(ind);

ind = ((c(1)*xf2 + c(2) - yf2)<0);
xf2u = xf2(ind);
yf2u = yf2(ind);

xfl = [xf1l; xf2l];
yfl = [yf1l; yf2l];

xfu = [xf1u; xf2u];
yfu = [yf1u; yf2u];

[xfl,ind] = sort(xfl);
yfl = yfl(ind); 

[xfu,ind] = sort(xfu);
yfu = yfu(ind); 

%figure(1);clf;plot(xfl,yfl,'-r',xfu,yfu,'-b'); axis equal; axis tight;
