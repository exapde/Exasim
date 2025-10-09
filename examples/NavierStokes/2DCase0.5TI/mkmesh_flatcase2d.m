function mesh = mkmesh_flatcase2d(porder, nx2, nx3, ny, dwall)
if nargin<2
    nx2 = 65;  % number of points on the flare line %orig 85
    nx3 = 9; %number of points on flat plate %orig 13
    ny = 71;   % number of points in the radial direction
    dwall = 0.0008; % control the element size closest to the wall        
end

s2 = 0.25;
s3 = 0.25; %flat plate addition
Lplate = 0.0298; %length of flat plate
%new cone segment
rnose = 0.0002;
thcone = 0;
rflare = 1;
Lcone = 0.3;
R = 0.06;
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
%yref = [0.8 0.4 0.15 0.06 0.015]*0.8; % locations to refine the mesh     
yref = [];

x2 = loginc(linspace(2,3,nx2), s2);
[p2,t2,yv] = lesmesh2d_rect(Href, dwall, ny, x2, yref);
 
%flat plate addition
x3 = loginc(linspace(1,2,nx3), s3);
[p3,t3] = lesmesh2d_rect(Href, dwall, ny, x3, yref);

[p,t] = connectmesh(p3',t3',p2',t2',1e-10);

elemtype = 1;
nodetype = 0;
bndexpr = {'all(p(:,2)<min(p0(:,2))+1e-6)','all(p(:,1)>max(p0(:,1))-1e-6)', ...
           'all(p(:,2)>max(p0(:,2))-1e-6)','true'};    
mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);
pnew = p;
ind3 = p(:,1)>1-1e-8 & p(:,1) <= 2+1e-8;
ind2 = p(:,1)>2+1e-8;          
 
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

ind3 = dgx > 1-1e-8 & dgx <= 2+1e-8;
ind2 = dgx > 2+1e-8;

%plate
dgx(ind3) = rnose + (dgx(ind3)-1)*Lplate;
dgy(ind3) = rnose + R * dgy(ind3);
%flare
x = thflare*(dgx(ind2)-2);
y = R*dgy(ind2);
dgx(ind2) = xc + (rflare-y).*sin(thcone + x);
dgy(ind2) = yc - (rflare-y).*cos(thcone + x);
  
%shift section
xshift = 0.02;
yshift = 0.0;

pnew(:,1) = pnew(:,1) + xshift;
pnew(:,2) = pnew(:,2) + yshift;

dgx = dgx + xshift;
dgy = dgy + yshift;

ind = pnew(:,2)==max(pnew(:,2));
p1 = pnew(ind,:);
ind = pnew(:,1)==max(pnew(:,1));
p2 = pnew(ind,:);
% [inflow, outflow, freestream, thermal wall]
bndexpr = {'all(p(:,1)<1e-6)', ...
           ['all((p(:,2)-' num2str(p1(2),12) ')/(' num2str(p1(2)-p2(2),12) ')' ' - (p(:,1)-' num2str(p1(1),12) ')/(' num2str(p1(1)-p2(1),12) ')>-1e-5)'],...            
           'all(p(:,2)>0.09-1e-6)',...                 
           'true'};            
mesh = mkmesh(pnew,t,porder,bndexpr,elemtype,nodetype);
mesh.dgnodes(:,1,:) = dgx;
mesh.dgnodes(:,2,:) = dgy;
 
mesh.boundaryexpr = {
  @(p) (p(1,:)< xshift + rnose + 1e-6),      ... %  Inflow 
  @(p) (p(2,:)< p2(2) + 1e-6), ... %  Wall  
  @(p) (p(1,:)> p1(1) - 1e-6),  ...%  Outflow
  @(p) (p(1,:) > -1e10)         %  Freestream
};

mesh.periodicboundary = [];
mesh.periodicexpr = {};

mesh.p = mesh.p';
mesh.t = mesh.t';
mesh.f = facenumbering(mesh.p,mesh.t,mesh.elemtype,mesh.boundaryexpr,mesh.periodicexpr);
mesh.xpe = mesh.plocal; 
mesh.telem = mesh.tlocal;

figure(1);clf;boundaryplot(mesh,1); hold on; %inlet
boundaryplot(mesh,2); %wall
boundaryplot(mesh,3); %outflow
boundaryplot(mesh,4); %freestream 

mesh.xpe = mesh.plocal;
mesh.telem = mesh.tlocal;

figure(2); clf; meshplot(mesh);
