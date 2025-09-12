function [p,t,dgnodes] = mkmesh_circincirc_Ma17b(porder,m,n,r1,r2,alpha)
% This routine is inspired from mkmesh_circincirc_half, but adapted to the
% Ma = 17.605 case presented in Barter Thesis.

elemtype = 1;
[p,t] = squaremesh(m-1,n-1,1,elemtype);
p=p'; 

%%%%%%%%%%% Defines 3 regions with different refining speeds %%%%%%%%%%%%%%
% Number of element in radial direction for each region
nn1 = ceil(5.5*n/10); nn2 = ceil(3.05*n/10); nn3 = n-nn1-nn2;
% Position of transition between refinement areas
x12 = 0.1; x23 = 0.375;
ln1 = loginc(linspace(0.00,x12,nn1+1),alpha);
ln2 = loginc(linspace(x12, x23,nn2+1),1);
ln3 = loginc(linspace(x23,1.00,nn3),2);
ln  = [ln1(1:end-1) ln2(1:end-1) ln3];
ln  = reshape(ones(m,1)*ln,[m*n,1]);

% Assign mesh point positions
p(:,2) = ln;
%p(:,2) = loginc(p(:,2),alpha); % OLD VERSION
%plot(ln,ln,'.'); % A SUPPRIMER
ind = p(:,1)<=0.5;
p(ind,1) = logdec(p(ind,1),1);
ind = p(:,1)>=0.5;
p(ind,1) = loginc(p(ind,1),1);

dgnodes = mkdgnodes(p',t,porder);

pnew = p;
pnew(:,1) = -(r1+(r2-r1)*p(:,2)).*sin(pi*p(:,1));
pnew(:,2) = -(r1+(r2-r1)*p(:,2)).*cos(pi*p(:,1));
[p,t] = fixmesh(pnew,t');
p = p';
t = t';

pnew = zeros(size(dgnodes));
pnew(:,1,:) = -(r1+(r2-r1)*dgnodes(:,2,:)).*sin(pi*dgnodes(:,1,:));
pnew(:,2,:) = -(r1+(r2-r1)*dgnodes(:,2,:)).*cos(pi*dgnodes(:,1,:));
dgnodes = pnew;


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
