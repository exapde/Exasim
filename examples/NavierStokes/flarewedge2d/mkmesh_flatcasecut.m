function mesh = mkmesh_flatcasecut(porder, nx2, ny, dwall)
 
if nargin<2
    nx2 = 75;  % number of points on the flare line %orig 85
    nx3 = 9; %number of points on flat plate %orig 13
    ny = 71;   % number of points in the radial direction
    dwall = 0.0002; % control the element size closest to the wall    
end
s2 = 0.25;%0.25;
s3 = 0.25; %flat plate addition
 
Lplate = 0.03; %length of flat plate (chop off by .02)
 
%new cone segment
rnose = 0.0002;
thcone = 0;
rflare = 1;
Lcone = 0.30;
R = 0.09;
 
%previous cone
% rnose = 0.00016;
% thcone = 1.5*pi/180;
% rflare = 3;
% Lcone = 0.56;
% R = 0.09;
% delta = 2e-5;
% xslip = 0.55;
 
xm = 0*(1-sin(thcone)) + Lplate; % move the center of the flared cone by Lplate
ym = rnose*cos(thcone);
xc = xm - rflare*sin(thcone);
yc = ym + rflare*cos(thcone);
xe = Lcone+Lplate;
ye = yc-sqrt(rflare^2 - (Lcone-xc)^2);
s = sqrt((xe-xm)^2+(ye-ym)^2);
thflare = acos(1 - s^2/(2*rflare^2));
 
Href = 1;  % reference height
yref = [0.5 0.25 0.1 0.04 0.015 0.004]*0.8; % locations to refine the mesh     
 
%x2 = loginc(linspace(2,3,nx2), s2);
%[p2,t2,yv] = lesmesh2d_rect(Href, dwall, ny, x2, yref);

% x1 = loginc(linspace(0,1,nx1),s1);
% [p1,t1] = quadgrid(x1,yv);
% alfa1 = asin(delta/rnose)/(pi/2-thcone);
% alfa2 = asin(delta/(rnose+R))/(pi/2-thcone);
% x1 = alfa1 + (alfa2-alfa1)*p1(:,2);
%alfa1 = asin(delta/rnose)/(pi/2-thcone);
% x1 = alfa1 + (alfa2-alfa1)*p1(:,2);
% alfa = asin(delta./(rnose+R*p1(:,2)))/(pi/2-thcone);
% p1(:,1) = alfa + (1-alfa).*p1(:,1);
 
%flat plate addition 
x3 = loginc(linspace(0,1,nx3), s3);
[p3,t3,yv] = lesmesh2d_rect(Href, dwall, ny, x3, yref);
% flare
x2 = loginc(linspace(1,2,nx2), s2);
[p2,t2, ~] = lesmesh2d_rect(Href, dwall, ny, x2, yref);
 
% connect plate -> flare
[p,t] = connectmesh(p3',t3',p2',t2',1e-10);
 
elemtype = 1;
nodetype = 0;
bndexpr = {'all(p(:,2)<min(p0(:,2))+1e-6)','all(p(:,1)>max(p0(:,1))-1e-6)', ...
           'all(p(:,2)>max(p0(:,2))-1e-6)','true'};     
mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);
 
pnew = p;
 
%ind1 = p(:,1)<=1+1e-8;            
ind3 = p(:,1) <= 1+1e-8;
ind2 = p(:,1) >  1+1e-8;
 
%need to offset by Lplate here?
%x = (pi/2-thcone)*p(ind1,1);
%y = R*p(ind1,2);
 
%nose
%pnew(ind1,1) = rnose - (rnose+y).*sqrt(1-(sin(x)).^2);
%pnew(ind1,2) = (rnose+y).*sin(x);
 
%plate
pnew(ind3,1) = (p(ind3,1))*Lplate;
pnew(ind3,2) = rnose + R*p(ind3,2);
 
%cone
x = thflare*(p(ind2,1)-1);
y = R*p(ind2,2);
pnew(ind2,1) = xc + (rflare-y).*sin(thcone + x);
pnew(ind2,2) = yc - (rflare-y).*cos(thcone + x);
 
%mesh.p = pnew;
%also need Lplate offsets here I think
dgx = mesh.dgnodes(:,1,:);
dgy = mesh.dgnodes(:,2,:);
 
ind3 = dgx <= 1+1e-8;   % plate
ind2 = dgx >  1+1e-8;   % flare

 
%nose
%x = (pi/2-thcone)*dgx(ind1);
%y = R*dgy(ind1);
%dgx(ind1) = rnose - (rnose+y).*sqrt(1-(sin(x)).^2);
%dgy(ind1) = (rnose+y).*sin(x);
 
%plate
dgx(ind3) =(dgx(ind3))*Lplate;
dgy(ind3) = rnose + R * dgy(ind3);
 
%flare
x = thflare*(dgx(ind2)-1);
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
  @(p) (p(1,:)< 1e-6),      ... %  Inflow 
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

