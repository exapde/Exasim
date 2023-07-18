function [x, y] = foilmapping(X, Y, func, param)
% (X, Y) : coordinates on the reference domain
% H      : the maximum height from the parametric curve
% func   : function that defines the paramatric curve and its derivatives
% param  : parameter inputs to the function func
% (x, y) : coordinates on the physical domain

% Example 1:
%   [p,t]= squaremesh(11,11,1,1);
%   [x, y] = normalmapping(p(:,1), p(:,2), 2, @halfcircle, []);
%   q = [x y];
%   figure(2);clf;simpplot(q,t); axis on; axis equal; axis tight;

% Example 2:
%   [p,t]= squaremesh(11,11,1,1); p(:,1) = p(:,1) - 1;
%   [x, y] = normalmapping(p(:,1), p(:,2), 5, @naca12, [1.0089304129]);
%   q = [x y];
%   figure(2);clf;simpplot(q,t); axis on; axis equal; axis tight;

H = param(1);

% parametric curve xf(X) and yf(X)
[xf, yf, dxf, dyf] = func(X, param(2:end));

% normal vector to the parametric curve
ndf = sqrt(dxf.^2 + dyf.^2);
nx =  -dyf./ndf;
ny =   dxf./ndf;

% nx = -dyf./ndf;
% ny =  dxf./ndf;

% fix the normal vector when the derivatives are unbounded
ii = find(abs(dyf(:))==Inf);
nx(ii) = -sign(dyf(ii));
ny(ii) = 0;

ii = find(abs(dxf(:))==Inf);
ny(ii) = -sign(dxf(ii));
nx(ii) = 0;

% distance between (xf, yf) and (x, y)
d = Y*H;

% coordinates on the physical domain
x = xf + d.*nx;
y = yf + d.*ny;

% figure(1); clf; 
% plot(xf, yf, 'o'); 
% hold on;
% plot(x, y, '*'); 
% axis equal; axis tight;
