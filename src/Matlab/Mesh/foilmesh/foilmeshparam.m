function [mesh, meshref] = foilmeshparam(porder, elemtype, nodetype, nxw, nxur, nxu, nxl, nr, sps, spr)
%CMESHPARAM  Creates mesh in parametric space for airfoil c-type grids
%
%        _________________________________________
%       |      |      |      |      |      |      |
%       |      |      |      |      |      |      |
%    nr |      |      |      |      |      |      |
%       |      |      |      |      |      |      |
%       |______|______|______|______|______|______|
%         nxw    nxur   nxuf   nxlf   nxlr   nxw
%
%   nxw  : number of subdivison in the wake
%   nxur : number of subdivision in the upper foil (rear)
%   nxuf : number of subdivision in the upper foil (front)
%   nxlf : number of subdivisions in the lower foil (front)
%   nxlr : number of subdivisions in the lower foil (rear)
%   nr   : number of subdivisions in the radial direction
%   sps(id) : streamwise size control
%     sps(1)  - ratio between the first and last elements in the wake 
%     sps(2)  - ratio between mid-chord and trailing edge element (lower)
%   spr(id) : radial size control
%     spr(1) - ratio between far-field and wake element (lower) 
%     spr(2) - ratio between far-field and trailing edge (lower)

% m = 
% [x1, y1] = ndgrid((0:m-1)/(m-1),(0:n-1)/(n-1));

% porder = 3;
% elemtype = 1;
% nodetype = 1;
H = 4;
L = 1.0089304129;
W = 4;

% nr = 5;
% nxw = 5;
% nxu = 5;
% nxl = 5;
% 
% sps = [3 8 3 8 1.5 1.5];
% spr = 0.01;
% spw = [3 3 3 3];

nr = 21;
nxw = 15;
nxu = 25;
nxl = 25;
sps = [3 8 3 8 3 3];
spr = 2;
spw = [3 3 3 3];

parity = 0;

% upper foil surface
mesh = mkmesh_square(nxu,nr,porder,parity,1,1,elemtype,nodetype);
pfu = mesh.p;
tfu = mesh.t;
dfu = mesh.dgnodes;
pfu(:,1) = pfu(:,1) - 1;
dfu(:,1,:) = dfu(:,1,:) - 1;

X = pfu(:,1); Y = pfu(:,2);
[x, y] = upperfoil(X, Y, sps(1), sps(2), spr(1), spw(1));
[x, y] = foilmapping(x, y, H, @naca12upper, L);
mfu = [x(:) y(:)];

X = dfu(:,1,:); Y = dfu(:,2,:);
[x, y] = upperfoil(X(:), Y(:), sps(1), sps(2), spr(1), spw(1));
[x, y] = foilmapping(x, y, H, @naca12upper, L);
nfu = 1.0*dfu;
nfu(:,1,:) = reshape(x, size(X));
nfu(:,2,:) = reshape(y, size(X));

% lower foil surface
mesh = mkmesh_square(nxl,nr,porder,parity,1,1,elemtype,nodetype);
pfl = mesh.p;
tfl = mesh.t;
dfl = mesh.dgnodes;

X = pfl(:,1); Y = pfl(:,2);
[x, y] = lowerfoil(X, Y, sps(3), sps(4), spr(1), spw(2));
[x, y] = foilmapping(x, y, H, @naca12lower, L);
mfl = [x(:) y(:)];

X = dfl(:,1,:); Y = dfl(:,2,:);
[x, y] = lowerfoil(X(:), Y(:), sps(3), sps(4), spr(1), spw(2));
[x, y] = foilmapping(x, y, H, @naca12lower, L);
nfl = 1.0*dfl;
nfl(:,1,:) = reshape(x, size(X));
nfl(:,2,:) = reshape(y, size(X));

% upper wake
mesh = mkmesh_square(nxw,nr,porder,parity,1,1,elemtype,nodetype);
pwu = mesh.p;
twu = mesh.t;
dwu = mesh.dgnodes;
pwu(:,1) = pwu(:,1) - 2;
dwu(:,1,:) = dwu(:,1,:) - 2;

X = pwu(:,1); Y = pwu(:,2);
[x, y] = upperwake(X, Y, sps(5), spr(1), spw(3));
q1 = mfu(1,:);
q2 = q1; q2(1) = q2(1)+W;
q3 = q2; q3(2) = q3(2)+H;
q4 = mfu(end-nxu+1,:);
mwu = mapp([x y], [q2; q1; q3; q4]);

X = dwu(:,1,:); Y = dwu(:,2,:);
[x, y] = upperwake(X(:), Y(:), sps(5), spr(1), spw(3));
nwu = mapp([x y], [q2; q1; q3; q4]);
nwu = permute(reshape(nwu, [size(dwu,1) size(dwu,3) 2]), [1 3 2]);

% lower wake
mesh = mkmesh_square(nxw,nr,porder,parity,1,1,elemtype,nodetype);
pwl = mesh.p;
twl = mesh.t;
dwl = mesh.dgnodes;
pwl(:,1) = pwl(:,1) + 1;
dwl(:,1,:) = dwl(:,1,:) + 1;

X = pwl(:,1); Y = pwl(:,2);
[x, y] = lowerwake(X, Y, sps(6), spr(1), spw(4));
q3 = q2; q3(2) = q3(2)-H;
q4 = mfl(end,:);
mwl = mapp([x y], [q1; q2; q4; q3]);

X = dwl(:,1,:); Y = dwl(:,2,:);
[x, y] = lowerwake(X(:), Y(:), sps(6), spr(1), spw(4));
nwl = mapp([x y], [q1; q2; q4; q3]);
nwl = permute(reshape(nwl, [size(dwl,1) size(dwl,3) 2]), [1 3 2]);

% connect subgrids to make the reference mesh 
[pref,tref] = connectmesh(pwu,twu,pfu,tfu,1e-6);   
[pref,tref] = connectmesh(pref,tref,pfl,tfl,1e-6); 
[pref,tref] = connectmesh(pref,tref,pwl,twl,1e-6);
[pref,tref] = fixmesh(pref,tref);
dgref = cat(3,dwu,dfu);
dgref = cat(3,dgref,dfl);
dgref = cat(3,dgref,dwl);

bndexpr = {'all(p(:,2)<1e-3)','all(p(:,1)>max(p0(:,1))-1e-3)', ...
           'all(p(:,2)>max(p0(:,2))-1e-3)', 'true'};     
meshref = mkmesh(pref,tref,porder,bndexpr,elemtype,nodetype);
meshref.dgnodes = dgref;

% connect subgrids to make the physical mesh 
[p,t] = connectmesh(mwu,twu,mfu,tfu,1e-6); 
[p,t] = connectmesh(p,t,mfl,tfl,1e-6); 
[p,t] = connectmesh(p,t,mwl,twl,1e-6);
[p,t] = fixmesh(p,t);
dgphy = cat(3,nwu,nfu);
dgphy = cat(3,dgphy,nfl);
dgphy = cat(3,dgphy,nwl);

bndexpr = {'all(sqrt(sum(p.^2,2))<2)','all(sqrt(sum(p.^2,2))>2)'};  
mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);
% mesh.dgnodes = dgphy;

% figure(1); clf; hold on;
% meshplot(meshref,[0 0 1]);
% axis equal; axis tight; axis on;
% 
% figure(2); clf; hold on;
% meshplot(mesh,[0 0 1]);
% axis equal; axis tight; axis on;

end


function [x, y] = upperfoil(X, Y, sps1, sps2, spr, c)

% warp X-coordinates to bend them toward the leading edge
a = min(X(:)); b = max(X(:)); 
x = a + (b-a)*(1-exp(-c*Y.*(X-a)/(b-a)))./(1-exp(-c*Y));
ii = abs(Y)<1e-8;
x(ii) = X(ii);

% warp Y-coordinates to bend them toward the leading edge
a = min(Y(:)); b = max(Y(:));
y = a + (b-a)*(exp(c*(1+X).*(Y-a)/(b-a))-1)./(exp(c*(1+X))-1);
ii = abs(1+X)<1e-8;
y(ii) = Y(ii);

% scale X-coordinates to make the grid fine at both the leading edge and the
% trailing edge
x = loginc(x, sps1);
x = logdec(x, sps2);

% scale Y-coordinates to make the grid fine near the foil surface
y = loginc(y, spr);

end

function [x, y] = lowerfoil(X, Y, sps1, sps2, spr, c)

% warp X-coordinates to bend them toward the leading edge
a = min(X(:)); b = max(X(:)); 
x = a + (b-a)*(exp(c*Y.*(X-a)/(b-a))-1)./(exp(c*Y)-1);
ii = abs(Y)<1e-8;
x(ii) = X(ii);

% warp Y-coordinates to bend them toward the leading edge
a = min(Y(:)); b = max(Y(:)); 
y = a + (b-a)*(exp(c*(1-X).*(Y-a)/(b-a))-1)./(exp(c*(1-X))-1);
ii = abs(1-X)<1e-8;
y(ii) = Y(ii);

% scale X-coordinates to make the grid fine at both the leading edge and the
% trailing edge
x = logdec(x, sps1);
x = loginc(x, sps2);

% scale Y-coordinates to make the grid fine near the foil surface
y = loginc(y, spr);

end

function [x, y] = upperwake(X, Y, sps, spr, c)

% warp X-coordinates to bend them toward the wake boundary
a = min(X(:)); b = max(X(:));
x = a + (b-a)*(exp(c*Y.*(X-a)/(b-a))-1)./(exp(c*Y)-1);
ii = abs(Y)<1e-8;
x(ii) = X(ii,1);

% scale X-coordinates to make the grid fine at the trailing edge
x = x + 2;
x = logdec(x, sps);

% scale Y-coordinates to make the grid fine along the wake centerline
y = loginc(Y, spr);

end

function [x, y] = lowerwake(X, Y, sps, spr, c)

% warp X-coordinates to bend them toward the wake boundary
a = min(X(:)); b = max(X(:));
x = a + (b-a)*(1-exp(-c*Y.*(X-a)/(b-a)))./(1-exp(-c*Y));
ii = abs(Y)<1e-8;
x(ii) = X(ii,1);

% scale X-coordinates to make the grid fine at the trailing edge
x = x - 1;
x = loginc(x, sps);

% scale Y-coordinates to make the grid fine along the wake centerline
y = loginc(Y, spr);

end

% figure(7); clf; hold on;
% simpplot(mfu, tfu); 
% simpplot(mfl, tfl); 
% simpplot(mwu, twu);
% simpplot(mwl, twl);
% axis equal; axis tight; axis on;

% figure(1); clf; hold on;
% simpplot(pfu, tfu);
% simpplot(pfl, tfl); 
% axis equal; axis tight; axis on;
% 
% figure(2); clf; hold on;
% simpplot(mfu, tfu); 
% simpplot(mfl, tfl); 
% axis equal; axis tight; axis on;

% a = min(X(:)); b = max(X(:)); c = 3;
% x = a + (b-a)*(1-exp(-c*Y.*(X-a)/(b-a)))./(1-exp(-c*Y));
% ii = abs(Y)<1e-8;
% x(ii) = pfu(ii,1);
% 
% a = min(Y(:)); b = max(Y(:)); c = 2;
% y = a + (b-a)*(exp(c*(1+X).*(Y-a)/(b-a))-1)./(exp(c*(1+X))-1);
% ii = abs(1+X)<1e-8;
% y(ii) = pfu(ii,2);
% 
% x = loginc(x, sps(2));
% x = logdec(x, sps(3));
% y = loginc(y, spr(1));



% a = min(X(:)); b = max(X(:)); c = 3;
% x = a + (b-a)*(exp(c*Y.*(X-a)/(b-a))-1)./(exp(c*Y)-1);
% ii = abs(Y)<1e-8;
% x(ii) = pfl(ii,1);
% 
% a = min(Y(:)); b = max(Y(:)); c = 2;
% y = a + (b-a)*(exp(c*(1-X).*(Y-a)/(b-a))-1)./(exp(c*(1-X))-1);
% ii = abs(1-X)<1e-8;
% y(ii) = pfl(ii,2);
% 
% x = logdec(x, sps(2));
% x = loginc(x, sps(3));
% y = loginc(y, spr(1));


% a = min(X(:)); b = max(X(:)); c = 3;
% x = a + (b-a)*(exp(c*Y.*(X-a)/(b-a))-1)./(exp(c*Y)-1);
% ii = abs(Y)<1e-8;
% x(ii) = X(ii,1);
% 
% x = x + 2;
% x = logdec(x, spr(1));
% y = Y;
% y = loginc(y, spr(1));


% a = min(X(:)); b = max(X(:)); c = 3;
% x = a + (b-a)*(1-exp(-c*Y.*(X-a)/(b-a)))./(1-exp(-c*Y));
% ii = abs(Y)<1e-8;
% x(ii) = X(ii,1);
% 
% x = x - 1;
% x = loginc(x, spr(1));
% y = Y;
% y = loginc(y, spr(1));
