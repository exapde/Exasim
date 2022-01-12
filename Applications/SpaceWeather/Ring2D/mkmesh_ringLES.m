function [p,t,dgnodes] = mkmesh_ringLES(porder,r1,r2)

dlay = 1;
dwall = 8e-3; 
nx = 121;
ny = 16; 
xv = linspace(0, 1, nx); 
yref = [3e-2 12e-2 3e-1 6e-1];

[p,t,yv] = lesmesh2d_rect(dlay, dwall, ny, xv, yref);
p = p';

dgnodes = mkdgnodes(p',t,porder);

pnew = p;
pnew(:,1) = -(r1+(r2-r1)*p(:,2)).*sin(2*pi*p(:,1));
pnew(:,2) = -(r1+(r2-r1)*p(:,2)).*cos(2*pi*p(:,1));
[p,t] = fixmesh(pnew,t');
p = p';
t = t';

pnew = zeros(size(dgnodes));
pnew(:,1,:) = -(r1+(r2-r1)*dgnodes(:,2,:)).*sin(2*pi*dgnodes(:,1,:));
pnew(:,2,:) = -(r1+(r2-r1)*dgnodes(:,2,:)).*cos(2*pi*dgnodes(:,1,:));
dgnodes = pnew;

% simpplot(p,t); axis equal; axis tight;

function dgnodes = mkdgnodes(p,t,porder)
%CREATEDGNODES Computes the Coordinates of the DG nodes.
%   DGNODES=CREATENODES(MESH,FD,FPARAMS)
%
%      MESH:      Mesh Data Structure
%      FD:        Distance Function d(x,y)
%      FPARAMS:   Additional parameters passed to FD
%      DGNODES:   Triangle indices (NPL,2,NT). The nodes on 
%                 the curved boundaries are projected to the
%                 true boundary using the distance function FD
%

% npv : number of nodes per volume element
% nfv : number of faces per volume element
% npf : number of nodes per face element

% if porder>4
%     error("app.porder must be less than or equal to 4.");
% end

[nve,ne]=size(t);
nd=size(p,1);

elemtype = 0;
if (nd==2) && (nve==4)
    elemtype=1;    
end
if (nd==3) && (nve==8)
    elemtype=1;    
end

plocal = masternodes(porder,nd,elemtype);

npl=size(plocal,1);
if nd==1
    xi  = plocal(:,1);
    philocal(:,1) = 1 - xi;
    philocal(:,2) = xi;
elseif nd==2 && nve==3 % tri
    xi  = plocal(:,1);
    eta = plocal(:,2);    
    philocal(:,1) = 1 - xi - eta;
    philocal(:,2) = xi;
    philocal(:,3) = eta;
elseif nd==2 && nve==4 % quad
    xi  = plocal(:,1);
    eta = plocal(:,2);
    philocal(:,1) = (1-xi).*(1-eta);
    philocal(:,2) = xi.*(1-eta);
    philocal(:,3) = xi.*eta;
    philocal(:,4) = (1-xi).*eta;
elseif nd==3 && nve==4 % tet
    xi   = plocal(:,1);
    eta  = plocal(:,2);
    zeta = plocal(:,3);
    philocal(:,1) = 1 - xi - eta - zeta;
    philocal(:,2) = xi;
    philocal(:,3) = eta;
    philocal(:,4) = zeta;
elseif nd==3 && nve==8 % hex
    xi   = plocal(:,1);
    eta  = plocal(:,2);
    zeta = plocal(:,3);
    philocal(:,1) = (1-xi).*(1-eta).*(1-zeta);
    philocal(:,2) = xi.*(1-eta).*(1-zeta);
    philocal(:,3) = xi.*eta.*(1-zeta);
    philocal(:,4) = (1-xi).*eta.*(1-zeta);    
    philocal(:,5) = (1-xi).*(1-eta).*(zeta);
    philocal(:,6) = xi.*(1-eta).*(zeta);
    philocal(:,7) = xi.*eta.*(zeta);
    philocal(:,8) = (1-xi).*eta.*(zeta);        
end
    
% Allocate nodes
dgnodes=zeros(npl,nd,ne);
for dim=1:nd
  for node=1:nve
    dp=reshape(philocal(:,node),[npl 1])*reshape(p(dim,t(node,:)),[1 ne]);
    dgnodes(:,dim,:)=dgnodes(:,dim,:)+reshape(dp,[npl 1 ne]);
  end
end
