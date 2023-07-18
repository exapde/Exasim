function mesh = foilmesh2(porder)
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

%porder = 5;
elemtype = 1;
nodetype = 1;
H = 4;
L = 1.0089304129;
W = 3;

% nr = 5;
% nxw = 5;
% nxu = 5;
% nxl = 5;
% 
% sps = [3 8 3 8 1.5 1.5];
% spr = 0.01;
% spw = [3 3 3 3];

nr = 15;
nxw = 11;
nxu = 21;
nxl = 21;
sps = [3 8 3 8 3 3];
spr = 4;
spw = [3 3 3 3];

parity = 0;
ufunc = @naca12upper;
uparam = [H L];
lfunc = @naca12lower;
lparam = [H L];
usps = [sps(1) sps(2) spw(1)];
lsps = [sps(3) sps(4) spw(2)];

% foil mesh
mesh = mkmesh_rect(nxu+nxl-1,nr,porder,parity,[-1 1 0 1],elemtype,nodetype);
tf = mesh.t;
pf = mesh.p;
df = mesh.dgnodes;

[x, y] = xyfoil(mesh.p(:,1), mesh.p(:,2), ufunc, uparam, usps, lfunc, lparam, lsps, spr);
mf = [x(:) y(:)];

X = mesh.dgnodes(:,1,:); Y = mesh.dgnodes(:,2,:);
[x, y] = xyfoil(X(:), Y(:), ufunc, uparam, usps, lfunc, lparam, lsps, spr);
nf = 1.0*mesh.dgnodes;
nf(:,1,:) = reshape(x, size(X));
nf(:,2,:) = reshape(y, size(X));

mesh = mkmesh(mf, tf, porder,{'true'},elemtype,nodetype);
figure(1); clf; hold on;
meshplot(mesh,[0 0 1]);
axis equal; axis tight; axis on;
 
% upper wake mesh
mesh = mkmesh_square(nxw,nr,porder,parity,1,1,elemtype,nodetype);
pwu = mesh.p;
twu = mesh.t;
dwu = mesh.dgnodes;
pwu(:,1) = pwu(:,1) + 1;
dwu(:,1,:) = dwu(:,1,:) + 1;

[x, y] = upperwake(pwu(:,1), pwu(:,2), sps(5), spr, spw(3));
q1 = mf(nxu+nxl-1,:);
q2 = q1; q2(1) = q2(1)+W;
q3 = q2; q3(2) = q3(2)+H;
q4 = mf(end,:);
mwu = mapp([x y], [q1; q2; q4; q3]);

X = dwu(:,1,:); Y = dwu(:,2,:);
[x, y] = upperwake(X(:), Y(:), sps(5), spr, spw(3));
nwu = mapp([x y], [q1; q2; q4; q3]);
nwu = permute(reshape(nwu, [size(dwu,1) size(dwu,3) 2]), [1 3 2]);

mesh = mkmesh(mwu,twu,porder,{'true'},elemtype,nodetype);
meshplot(mesh,[0 0 1]);

% lower wake mesh
mesh = mkmesh_square(nxw,nr,porder,parity,1,1,elemtype,nodetype);
pwl = mesh.p;
twl = mesh.t;
dwl = mesh.dgnodes;
pwl(:,1) = pwl(:,1) - 2;
dwl(:,1,:) = dwl(:,1,:) - 2;

[x, y] = lowerwake(pwl(:,1), pwl(:,2), sps(6), spr, spw(4));
q3 = q2; q3(2) = q3(2)-H;
q4 = mf(end-nxu-nxl+2,:);
mwl = mapp([x y], [q2; q1; q3; q4]);

X = dwl(:,1,:); Y = dwl(:,2,:);
[x, y] = lowerwake(X(:), Y(:), sps(6), spr, spw(4));
nwl = mapp([x y], [q2; q1; q3; q4]);
nwl = permute(reshape(nwl, [size(dwl,1) size(dwl,3) 2]), [1 3 2]);

mesh = mkmesh(mwl,twl,porder,{'true'},elemtype,nodetype);
meshplot(mesh,[0 0 1]);

% connect subgrids to make the reference mesh 
[pref,tref] = connectmesh(pwl,twl,pf,tf,1e-4);   
[pref,tref] = connectmesh(pref,tref,pwu,twu,1e-4); 
%[pref,tref] = fixmesh(pref,tref);
dgref = cat(3,dwl,df);
dgref = cat(3,dgref,dwu);

bndexpr = {'all(p(:,2)<1e-3)','all(p(:,1)>max(p0(:,1))-1e-3)', ...
           'all(p(:,2)>max(p0(:,2))-1e-3)', 'true'};     
meshref = mkmesh(pref,tref,porder,bndexpr,elemtype,nodetype);
meshref.dgnodes = dgref;

% connect subgrids to make the physical mesh 
[p,t] = connectmesh(mwl,twl,mf,tf,1e-4); 
[p,t] = connectmesh(p,t,mwu,twu,1e-4); 
%[p,t] = fixmesh(p,t);
dgphy = cat(3,nwl,nf);
dgphy = cat(3,dgphy,nwu);

bndexpr = {'all(sqrt(sum(p.^2,2))<2)','all(sqrt(sum(p.^2,2))>2)'};  
mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);
mesh.dgnodes = dgphy;

figure(2); clf; hold on;
meshplot(meshref,0);
axis equal; axis tight; axis on;

figure(3); clf; hold on;
meshplot(mesh,0);
axis equal; axis tight; axis on;


% % upper foil surface
% mesh = mkmesh_square(nxu,nr,porder,parity,1,1,elemtype,nodetype);
% pfu = mesh.p;
% tfu = mesh.t;
% dfu = mesh.dgnodes;
% 
% X = pfu(:,1); Y = pfu(:,2);
% [x, y] = upperfoil(X, Y, @naca12upper, [H L], [sps(1) sps(2) spw(1)], spr);
% mfu = [x(:) y(:)];
% 
% X = dfu(:,1,:); Y = dfu(:,2,:);
% [x, y] = upperfoil(X(:), Y(:), @naca12upper, [H L], [sps(1) sps(2) spw(1)], spr);
% nfu = 1.0*dfu;
% nfu(:,1,:) = reshape(x, size(X));
% nfu(:,2,:) = reshape(y, size(X));
% 
% mesh = mkmesh(mfu,tfu,porder,{'true'},elemtype,nodetype);
% figure(1); clf; hold on;
% meshplot(mesh,[0 0 1]);
% axis equal; axis tight; axis on;
% 
% % % lower foil surface
% mesh = mkmesh_square(nxl,nr,porder,parity,1,1,elemtype,nodetype);
% pfl = mesh.p;
% tfl = mesh.t;
% dfl = mesh.dgnodes;
% pfl(:,1) = pfl(:,1) - 1;
% dfl(:,1,:) = dfl(:,1,:) - 1;
% 
% X = pfl(:,1); Y = pfl(:,2);
% [x, y] = lowerfoil(X, Y, @naca12lower, [H L], [sps(3), sps(4), spw(2)], spr);
% mfl = [x(:) y(:)];
% 
% X = dfl(:,1,:); Y = dfl(:,2,:);
% [x, y] = lowerfoil(X(:), Y(:), @naca12lower, [H L], [sps(3), sps(4), spw(2)], spr);
% nfl = 1.0*dfl;
% nfl(:,1,:) = reshape(x, size(X));
% nfl(:,2,:) = reshape(y, size(X));
% 
% mesh = mkmesh(mfl,tfl,porder,{'true'},elemtype,nodetype);
% %figure(1); clf; hold on;
% meshplot(mesh,[0 0 1]);
% 
% % upper wake
% mesh = mkmesh_square(nxw,nr,porder,parity,1,1,elemtype,nodetype);
% pwu = mesh.p;
% twu = mesh.t;
% dwu = mesh.dgnodes;
% pwu(:,1) = pwu(:,1) + 1;
% dwu(:,1,:) = dwu(:,1,:) + 1;
% 
% X = pwu(:,1); Y = pwu(:,2);
% [x, y] = upperwake(X, Y, sps(5), spr(1), spw(3));
% q1 = mfu(nxu,:);
% q2 = q1; q2(1) = q2(1)+W;
% q3 = q2; q3(2) = q3(2)+H;
% q4 = mfu(end,:);
% mwu = mapp([x y], [q1; q2; q4; q3]);
% 
% X = dwu(:,1,:); Y = dwu(:,2,:);
% [x, y] = upperwake(X(:), Y(:), sps(5), spr(1), spw(3));
% nwu = mapp([x y], [q2; q1; q3; q4]);
% nwu = permute(reshape(nwu, [size(dwu,1) size(dwu,3) 2]), [1 3 2]);
% 
% mesh = mkmesh(mwu,twu,porder,{'true'},elemtype,nodetype);
% meshplot(mesh,[0 0 1]);
% 
% % lower wake
% mesh = mkmesh_square(nxw,nr,porder,parity,1,1,elemtype,nodetype);
% pwl = mesh.p;
% twl = mesh.t;
% dwl = mesh.dgnodes;
% pwl(:,1) = pwl(:,1) + 1;
% dwl(:,1,:) = dwl(:,1,:) + 1;
% 
% X = pwl(:,1); Y = pwl(:,2);
% [x, y] = lowerwake(X, Y, sps(6), spr(1), spw(4));
% q3 = q2; q3(2) = q3(2)-H;
% q4 = mfl(end,:);
% mwl = mapp([x y], [q1; q2; q4; q3]);
% 
% X = dwl(:,1,:); Y = dwl(:,2,:);
% [x, y] = lowerwake(X(:), Y(:), sps(6), spr(1), spw(4));
% nwl = mapp([x y], [q1; q2; q4; q3]);
% nwl = permute(reshape(nwl, [size(dwl,1) size(dwl,3) 2]), [1 3 2]);
% 
% % connect subgrids to make the reference mesh 
% [pref,tref] = connectmesh(pwu,twu,pfu,tfu,1e-6);   
% [pref,tref] = connectmesh(pref,tref,pfl,tfl,1e-6); 
% [pref,tref] = connectmesh(pref,tref,pwl,twl,1e-6);
% [pref,tref] = fixmesh(pref,tref);
% dgref = cat(3,dwu,dfu);
% dgref = cat(3,dgref,dfl);
% dgref = cat(3,dgref,dwl);
% 
% bndexpr = {'all(p(:,2)<1e-3)','all(p(:,1)>max(p0(:,1))-1e-3)', ...
%            'all(p(:,2)>max(p0(:,2))-1e-3)', 'true'};     
% meshref = mkmesh(pref,tref,porder,bndexpr,elemtype,nodetype);
% meshref.dgnodes = dgref;
% 
% % connect subgrids to make the physical mesh 
% [p,t] = connectmesh(mwu,twu,mfu,tfu,1e-6); 
% [p,t] = connectmesh(p,t,mfl,tfl,1e-6); 
% [p,t] = connectmesh(p,t,mwl,twl,1e-6);
% [p,t] = fixmesh(p,t);
% dgphy = cat(3,nwu,nfu);
% dgphy = cat(3,dgphy,nfl);
% dgphy = cat(3,dgphy,nwl);
% 
% bndexpr = {'all(sqrt(sum(p.^2,2))<2)','all(sqrt(sum(p.^2,2))>2)'};  
% mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);
% mesh.dgnodes = dgphy;

% figure(1); clf; hold on;
% meshplot(meshref,[0 0 1]);
% axis equal; axis tight; axis on;
% 
% figure(2); clf; hold on;
% meshplot(mesh,[0 0 1]);
% axis equal; axis tight; axis on;

end

function [x, y] = xyfoil(X, Y, ufunc, uparam, usps, lfunc, lparam, lsps, spr)

iu = X >= 0;
il = X <= 0;

x = 0*X; y = 0*Y;
 
[x1, y1] = upperfoil(X(iu), Y(iu), ufunc, uparam, usps, spr);
x(iu) = x1; y(iu) = y1; 

[x1, y1] = lowerfoil(X(il), Y(il), lfunc, lparam, lsps, spr);
x(il) = x1; y(il) = y1; 

end

function [x, y] = upperfoil(X, Y, upperfunc, param, sps, spr)

sps1 = sps(1);
sps2 = sps(2);
c = sps(3);

% warp X-coordinates to bend them toward the leading edge
%a = min(X(:)); b = max(X(:)); 
a=0; b=1;
x = a + (b-a)*(exp(c*Y.*(X-a)/(b-a))-1)./(exp(c*Y)-1);
ii = abs(Y)<1e-8;
x(ii) = X(ii);
 
% warp Y-coordinates to bend them toward the leading edge
a = min(Y(:)); b = max(Y(:));
y = a + (b-a)*(exp(c*(1-X).*(Y-a)/(b-a))-1)./(exp(c*(1-X))-1);
ii = abs(1-X)<1e-8;
y(ii) = Y(ii);

% x = X; y = Y;

% scale X-coordinates to make the grid fine at both the leading edge and the
% trailing edge
% x = logdec(x, sps1);
% x = loginc(x, sps2);
a=0; b=1;
alpha = sps1;
x = a + (b-a)*(1-exp(-alpha*(x-a)/(b-a)))/(1-exp(-alpha));
alpha = sps2;
x = a + (b-a)*(exp(alpha*(x-a)/(b-a))-1)/(exp(alpha)-1);

% scale Y-coordinates to make the grid fine near the foil surface
y = loginc(y, spr);

[x, y] = foilmapping(x, y, upperfunc, param);

end

function [x, y] = lowerfoil(X, Y, lowerfunc, param, sps, spr)

sps1 = sps(1);
sps2 = sps(2);
c = sps(3);

% warp X-coordinates to bend them toward the leading edge
%a = min(X(:)); b = max(X(:)); 
a=-1; b=0;
x = a + (b-a)*(1-exp(-c*Y.*(X-a)/(b-a)))./(1-exp(-c*Y));
ii = abs(Y)<1e-8;
x(ii) = X(ii);

% warp Y-coordinates to bend them toward the leading edge
a = min(Y(:)); b = max(Y(:)); 
y = a + (b-a)*(exp(c*(X+1).*(Y-a)/(b-a))-1)./(exp(c*(X+1))-1);
ii = abs(X+1)<1e-8;
y(ii) = Y(ii);

% x = X; y = Y;

% scale X-coordinates to make the grid fine at both the leading edge and the
% trailing edge
% x = loginc(x, sps1);
% x = logdec(x, sps2);
a=-1; b=0;
alpha = sps1;
x = a + (b-a)*(exp(alpha*(x-a)/(b-a))-1)/(exp(alpha)-1);
alpha = sps2;
x = a + (b-a)*(1-exp(-alpha*(x-a)/(b-a)))/(1-exp(-alpha));

% scale Y-coordinates to make the grid fine near the foil surface
y = loginc(y, spr);

[x, y] = foilmapping(x, y, lowerfunc, param);

end

function [x, y] = upperwake(X, Y, sps, spr, c)

% warp X-coordinates to bend them toward the wake boundary
a = min(X(:)); b = max(X(:));
%x = a + (b-a)*(exp(c*Y.*(X-a)/(b-a))-1)./(exp(c*Y)-1);
x = a + (b-a)*(1-exp(-c*Y.*(X-a)/(b-a)))./(1-exp(-c*Y));
ii = abs(Y)<1e-8;
x(ii) = X(ii,1);

% x = X; y = Y;

% scale X-coordinates to make the grid fine at the trailing edge
x = x - 1;
x = loginc(x, sps);

% scale Y-coordinates to make the grid fine along the wake centerline
y = loginc(Y, spr);

end

function [x, y] = lowerwake(X, Y, sps, spr, c)

% warp X-coordinates to bend them toward the wake boundary
a = min(X(:)); b = max(X(:));
%x = a + (b-a)*(1-exp(-c*Y.*(X-a)/(b-a)))./(1-exp(-c*Y));
x = a + (b-a)*(exp(c*Y.*(X-a)/(b-a))-1)./(exp(c*Y)-1);
ii = abs(Y)<1e-8;
x(ii) = X(ii,1);

% x = X; y = Y;

% scale X-coordinates to make the grid fine at the trailing edge
x = x + 2;
x = logdec(x, sps);

% scale Y-coordinates to make the grid fine along the wake centerline
y = loginc(Y, spr);

end

