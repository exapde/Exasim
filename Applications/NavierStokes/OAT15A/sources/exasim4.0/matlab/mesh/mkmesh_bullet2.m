function mesh = mkmesh_bullet2(porder,elemtype,gridNum)

if nargin<1, porder=1;   end
if nargin<2, elemtype=1; end
if nargin<3, gridNum=1;  end

if gridNum==1
   n1=18*porder+1; n2=12*porder+1; n3=18*porder+1; 
   [x,y] = cmeshparam6(n1, n2, n2, n2, n2, n3, ...
                        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], ...
                        [10, 10, 10, 10, 10, 10, 10]*5);   
elseif gridNum==2
   n1=12*porder+1; n2=10*porder+1; n3=18*porder+1; 
   [x,y] = cmeshparam6(n1, n2, n2, n2, n2, n3, ...
                        [5, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1], ...
                        [10, 10, 10, 10, 10, 10, 10]*1e4);      
elseif gridNum==3
   n1=14*porder+1; n2=6*porder+1; n3=10*porder+1;    
   [x,y] = cmeshparam6(n1, n2, n2, n2, n2, n3, ...
                        [20/4, 20/4, 20/4, 20/4, 20/4, 1, 1, 1, 1, 1, 1], ...
                        [10, 10, 10, 10, 10, 10, 10]*5);   
end

[xu, yu, xl, yl] = surfacegeo;

% th = (pi:-pi/500:pi/2)';
% xt = (cos(th)+1);
% xu = xt(end:-1:1);
% xl = flipud(xu);   
% yu = sqrt(0.25 - (xu-0.5).^2);
% yl = -sqrt(0.25 - (xl-0.5).^2);

xf = [xu; xl(2:end)];
yf = [yu; yl(2:end)];
% xf(end) = xf(1);
% yf(end) = yf(1);

figure(1); clf;
plot(xf, yf, '-b');
axis equal;
axis tight;

[xm, ym] = cmeshmap(xf, yf, x, y, 6, 4);
% fix the wake gap
xm(1,1:n1) = xm(1,end:-1:end-(n1-1));
ym(1,1:n1) = ym(1,end:-1:end-(n1-1));

bndexpr={'sqrt((p(:,1)-.5).^2+p(:,2).^2)<2','true'};
mesh = cart2mesh(porder,xm,ym,[],bndexpr,elemtype);

function [xu, yu, xl, yl] = surfacegeo

L = 1;    % length of the vehicle
a = 0.2*L;
b = 0.2*L;
c = 0.2*L;
d = 0.4*L;
r = (b^2+d^2)/(2*b);
e = 2*a;
f = c/sqrt(1 - (e-a)^2/e^2);

x0 = 0;
x1 = a;
x2 = L-b;
x3 = L;

n  = 100;

xu1 = linspace(x3,x2,n)';
yu1 = sqrt(r^2 - (xu1-(L-r)).^2);
xu2 = linspace(x2,x1,n)';
yu2 = ((d-c)/(L-a-b))*xu2 + c - (d-c)*a/(L-a-b);
xu3 = linspace(x1,x0,n)';
yu3 = f*sqrt(1 - (xu3-e).^2/e^2);
xu  = [xu1; xu2(2:end); xu3(2:end)];
yu  = [yu1; yu2(2:end); yu3(2:end)];

xl1 = linspace(x0,x1,n)';
yl1 = -f*sqrt(1 - (xl1-e).^2/e^2);
xl2 = linspace(x1,x2,n)';
yl2 = -((d-c)/(L-a-b))*xl2 - c + (d-c)*a/(L-a-b);
xl3 = linspace(x2,x3,n)';
yl3 = -sqrt(r^2 - (xl3-(L-r)).^2);
xl  = [xl1; xl2(2:end); xl3(2:end)];
yl  = [yl1; yl2(2:end); yl3(2:end)];

% figure(1); clf;
% plot(xu, yu, '-b', xl, yl, '--r');
% axis equal;
% axis tight;
% pause
% yu = sqrt(0.25 - (xu-0.5).^2);
% yl = -sqrt(0.25 - (xl-0.5).^2);

