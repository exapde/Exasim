function mesh = mkmesh_flaredplate2d(porder, nx1, nx2, ny, dwall)
if nargin<2
    nx1 = 21;  % number of points on the nose line
    nx2 = 101;  % number of points on the flare line
    nx3 = 51; %number of points on flat plate
    ny = 101;   % number of points in the radial direction
    dwall = 0.0003; % control the element size closest to the wall   
end
s1 = 0.5;
s2 = 0.1;
s3 = 3.5; %flat plate addition
Lplate = 0.0498; %length of flat plate
%new cone segment
rnose = 0.0002;
thcone = 0;
rflare = 1;
Lcone = 0.35;
R = 0.1;
delta = 0;

xm = rnose*(1-sin(thcone)) + Lplate; % move the center of the flared cone by Lplate
ym = rnose*cos(thcone);
xc = xm - rflare*sin(thcone);
yc = ym + rflare*cos(thcone);
%xe = Lcone;
%ye = yc-sqrt(rflare^2 - (Lcone-xc)^2);
xe = Lplate + Lcone; %new
ye = yc - sqrt(rflare^2 - (xe - xc)^2);%new
s = sqrt((xe - xm)^2 + (ye - ym)^2);%new
thflare = acos(1 - s^2/(2*rflare^2));
Href = 1;  % reference height
yref = [];

x2 = loginc(linspace(2,3,nx2), s2);
[p2,t2,yv] = lesmesh2d_rect(Href, dwall, ny, x2, yref);
 
x1 = loginc(linspace(0,1,nx1),s1);
[p1,t1] = quadgrid(x1,yv);
alfa = asin(delta./(rnose+R*p1(:,2)))/(pi/2-thcone);
p1(:,1) = alfa + (1-alfa).*p1(:,1);
%flat plate addition
x3 = loginc(linspace(1,2,nx3), s3);
[p3,t3] = lesmesh2d_rect(Href, dwall, ny, x3, yref);
[p13, t13] = connectmesh(p1,t1,p3',t3',1e-10);
[p,t] = connectmesh(p13,t13,p2',t2',1e-10);
 
elemtype = 1;
nodetype = 0;
bndexpr = {'all(p(:,2)<min(p0(:,2))+1e-6)','all(p(:,1)>max(p0(:,1))-1e-6)', ...
           'all(p(:,2)>max(p0(:,2))-1e-6)','true'};    
mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);
pnew = p;
ind1 = p(:,1)<=1+1e-8;           
ind3 = p(:,1)>1-1e-8 & p(:,1) <= 2+1e-8;
ind2 = p(:,1)>2+1e-8;          
 
%need to offset by Lplate here?
x = (pi/2-thcone)*p(ind1,1);
y = R*p(ind1,2);
%nose
pnew(ind1,1) = rnose - (rnose+y).*sqrt(1-(sin(x)).^2);
pnew(ind1,2) = (rnose+y).*sin(x);
%plate
pnew(ind3,1) = rnose + (p(ind3,1)-1)*Lplate;
pnew(ind3,2) = rnose + R*p(ind3,2);
%cone
x = thflare*(p(ind2,1)-2);
y = R*p(ind2,2);
pnew(ind2,1) = xc + (rflare-y).*sin(thcone + x);
pnew(ind2,2) = yc - (rflare-y).*cos(thcone + x);

%mesh.p = pnew;
dgx = mesh.dgnodes(:,1,:);
dgy = mesh.dgnodes(:,2,:);
ind1 = dgx <= 1+1e-8;
ind3 = dgx > 1+1e-8 & dgx <= 2+1e-8;
ind2 = dgx > 2+1e-8;
%nose
x = (pi/2-thcone)*dgx(ind1);
y = R*dgy(ind1);
dgx(ind1) = rnose - (rnose+y).*sqrt(1-(sin(x)).^2);
dgy(ind1) = (rnose+y).*sin(x);
%plate
dgx(ind3) = rnose + (dgx(ind3)-1)*Lplate;
dgy(ind3) = rnose + R * dgy(ind3);
%flare
x = thflare*(dgx(ind2)-2);
y = R*dgy(ind2);
dgx(ind2) = xc + (rflare-y).*sin(thcone + x);
dgy(ind2) = yc - (rflare-y).*cos(thcone + x);
  
xO = 1e-6;   yO = delta;
xL = rnose; yL = rnose;% leading edge of flat plate
xC = rnose+Lplate;  yC = rnose;               % leading edge of flared] cone
%xC = xc;     yC = yc;           % leading edge of flared cone
xE = xe;     yE = ye;              % end point of flared cone
yT = max(pnew(:,2));  % top?right corner
idx = pnew(:,2)==yT;
xT = pnew(idx,1);
xF = .3498;
bndexpr = {'true'};
 % [xO yO]
% [xL yL]
% [xC yC]
% [xE yE]
% [xT yT]
 
mesh = mkmesh(pnew,t,porder,bndexpr,elemtype,nodetype);


mesh.dgnodes(:,1,:) = dgx; 
mesh.dgnodes(:,2,:) = dgy;
 
offset = .10625;
mesh.dgnodes(:,2,:) = mesh.dgnodes(:,2,:) + offset; 
mesh.p(:,2) = mesh.p(:,2) + offset;
yO = yO + offset;
yL = yL + offset;
yC = yC + offset;
yE = yE + offset;

mesh.boundaryexpr = {@(p) (p(2,:)<yO + 1e-6 & p(1,:)<xO + 1e-6),         ... % 1 Slip wall
  @(p) (p(1,:)<xL + 1e-6) &  sqrt(p(1,:).^2+(p(2,:)-offset).^2)>2*rnose+1e-6,                           ... % 2 Inflow (cap arc)
  @(p) (p(1,:)<xL + 1e-6) &  sqrt(p(1,:).^2+(p(2,:)-offset).^2)<2*rnose+1e-6,                           ... % 3 cap
  @(p) (p(2,:)<yL + 1e-6 & p(1,:)>xL - 1e-6 & p(1,:)<xC + 1e-6), ... % 4 Flat plate
  @(p) (p(2,:)<yE + 1e-6 & p(1,:)>xC - 1e-6 & p(1,:)<xF + 1e-6),         ... % 5 Flared cone (noslip)
  @(p) (p(2,:)<yE + 1e-6 & p(1,:)>xF-.1 - 1e-6 & p(1,:)<xE + 1e-6),         ... % 6 Flared cone (slip)
  @(p) (p(1,:)>xT - 1e-6),  ... % 7 Outflow
  @(p) (p(1,:) > -1e10)                                          ... % 8 Freestream
};
  
mesh.periodicboundary = [];
mesh.periodicexpr = {};
%mesh.boundarycondition = [1, 1, 2, 3, 1, 1];
 
mesh.p = mesh.p';
mesh.t = mesh.t';

mesh.f = facenumbering(mesh.p,mesh.t,mesh.elemtype,mesh.boundaryexpr,mesh.periodicexpr);
mesh.xpe = mesh.plocal; 
 
figure(1);clf;boundaryplot(mesh,1); hold on;
boundaryplot(mesh,2);boundaryplot(mesh,3);
boundaryplot(mesh,4);boundaryplot(mesh,5);
boundaryplot(mesh,6);boundaryplot(mesh,7);
boundaryplot(mesh,8);
plot(xO, yO, 'or', 'Linewidth', 2, 'MarkerSize', 8);
plot(xL, yL, 'or', 'Linewidth', 2, 'MarkerSize', 8);
plot(xC, yC, 'or', 'Linewidth', 2, 'MarkerSize', 8);
plot(xE, yE, 'or', 'Linewidth', 2, 'MarkerSize', 8);
plot(xT, yT, 'or', 'Linewidth', 2, 'MarkerSize', 8);
axis equal; axis tight;
 
mesh.xpe = mesh.plocal;
mesh.telem = mesh.tlocal;


