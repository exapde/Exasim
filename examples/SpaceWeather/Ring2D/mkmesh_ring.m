function [p,t,dgnodes] = mkmesh_ring(porder,m,n,r1,r2,alpha)
% This routine is inspired from mkmesh_circincirc_half, but adapted to the
% geometry of a ring.

elemtype = 1;
[p,t] = squaremesh(m-1,n-1,1,elemtype);
p=p'; 

ln = loginc(linspace(0.00,1,n),alpha);
ln  = reshape(ones(m,1)*ln,[m*n,1]);

% Assign mesh point positions
p(:,2) = ln;

ind = p(:,1)<=0.5;
p(ind,1) = logdec(p(ind,1),1);
ind = p(:,1)>=0.5;
p(ind,1) = loginc(p(ind,1),1);

dgnodes = mkdgnodes(p',t,porder);

pnew = p;
pnew(:,1) = -(r1+(r2-r1)*p(:,2)).*sin(2*pi*p(:,1));
pnew(:,2) = -(r1+(r2-r1)*p(:,2)).*cos(2*pi*p(:,1));
[p,t] = fixmeshCirc(pnew,t');
p = p';
t = t';

pnew = zeros(size(dgnodes));
pnew(:,1,:) = -(r1+(r2-r1)*dgnodes(:,2,:)).*sin(2*pi*dgnodes(:,1,:));
pnew(:,2,:) = -(r1+(r2-r1)*dgnodes(:,2,:)).*cos(2*pi*dgnodes(:,1,:));
dgnodes = pnew;


function [p,t]=fixmeshCirc(p,t)
%FIXMESH  Remove duplicated/unused nodes and fix element orientation.
%   [P,T]=FIXMESH(P,T)
 
% NOTE: This function works well for triangles, quadrangles and tetraedra.
% For hexaedra, this function can only recover elements whose bottom and up
% faces are numbered in a wrong (clockwise) way. TANGLED HEXES are not
% fixed with this routine !!!
 
disp('Fixing mesh...')
 
% Remove duplicated nodes:
snap=max(max(p,[],1)-min(p,[],1),[],2)*1024*eps;
[foo,ix,jx]=unique(round(p/snap)*snap,'rows');
if size(p,1) ~= length(ix); warning('Some vertices in mesh.p have been removed.'); end
p=p(ix,:);
t=jx(t);
if size(t,2) == 1, t = t'; end  % This lines ensures the function works for one element
 
% Remove nodes that are not contained in t:
[pix,ix,jx]=unique(t);
t=reshape(jx,size(t));
p=p(pix,:);
 
if (size(t,2) == 3 && size(p,2) == 2) || (size(t,2) == 4 && size(p,2) == 3)          % Simplices
    v = simpvol(p,t);
    flip=v<0;
    t(flip,[1,2])=t(flip,[2,1]);
elseif (size(t,2) == 4 && size(p,2) == 2)      % Quads
    D1 = p(t(:,3),:) - p(t(:,1),:);
    D2 = p(t(:,4),:) - p(t(:,2),:);
    flip=(D1(:,1).*D2(:,2) - D1(:,2).*D2(:,1))<0;
    t(flip,[1,2,3,4])=t(flip,[4,3,2,1]);
elseif (size(t,2) == 8 && size(p,2) == 3)      % Hex
    V12 = p(t(:,2),:) - p(t(:,1),:);
    V14 = p(t(:,4),:) - p(t(:,1),:);
    VAC = 1./4. * (p(t(:,5),:) + p(t(:,6),:) + p(t(:,7),:) + p(t(:,8),:)...
                 -(p(t(:,1),:) + p(t(:,2),:) + p(t(:,3),:) + p(t(:,4),:)));
    N1 = cross(V12,V14,2);
    flip = dot(N1,VAC,2)<0.;
    t(flip,[1,2,3,4])=t(flip,[4,3,2,1]);
    t(flip,[5,6,7,8])=t(flip,[8,7,6,5]);
    fprintf('%d Hexes have been reordered.\n',sum(flip))
else
    error('fixmesh not valid for this type of elements.');
end
 
if any(flip); warning('Some vertices in mesh.t have been reordered to meet code requirements.'); end



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
