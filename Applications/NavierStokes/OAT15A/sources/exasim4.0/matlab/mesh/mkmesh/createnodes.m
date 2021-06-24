function dgnodes=createnodes(mesh,fd,varargin)
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

if nargin < 2, fd=[]; end

p = mesh.p;
t = mesh.t;
plocal = mesh.plocal;

npl=size(plocal,1);
[nt,npv]=size(t);
nd=size(p,2);

if nd==1
    xi  = plocal(:,1);
    philocal(:,1) = 1 - xi;
    philocal(:,2) = xi;
elseif nd==2 && npv==3 % tri
    xi  = plocal(:,1);
    eta = plocal(:,2);    
    philocal(:,1) = 1 - xi - eta;
    philocal(:,2) = xi;
    philocal(:,3) = eta;
elseif nd==2 && npv==4 % quad
    xi  = plocal(:,1);
    eta = plocal(:,2);
    philocal(:,1) = (1-xi).*(1-eta);
    philocal(:,2) = xi.*(1-eta);
    philocal(:,3) = xi.*eta;
    philocal(:,4) = (1-xi).*eta;
elseif nd==3 && npv==4 % tet
    xi   = plocal(:,1);
    eta  = plocal(:,2);
    zeta = plocal(:,3);
    philocal(:,1) = 1 - xi - eta - zeta;
    philocal(:,2) = xi;
    philocal(:,3) = eta;
    philocal(:,4) = zeta;
elseif nd==3 && npv==8 % hex
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
dgnodes=zeros(npl,nd,nt);
for dim=1:nd
  for node=1:npv
    dp=philocal(:,node)*p(t(:,node),dim)';
    dgnodes(:,dim,:)=dgnodes(:,dim,:)+permute(dp,[1,3,2]);
  end
end

% Project nodes on the curved boundary
if ~isempty(fd) && mesh.porder>1
  tc=find(mesh.tcurved);
  if nd == 2
      for it=tc'
        p = dgnodes(:,:,it);
        deps=sqrt(eps)*max(max(p)-min(p));
        ed = find(mesh.f(abs(mesh.t2f(it,:)),end)<0);
        for id=ed'        
            e = find(philocal(:,id) < 1.e-6);        
            %d=feval(fd,p(e,:),varargin{:});        
            %dgradx=(feval(fd,[p(e,1)+deps,p(e,2)],varargin{:})-d)/deps;
            %dgrady=(feval(fd,[p(e,1),p(e,2)+deps],varargin{:})-d)/deps;
            d=feval(fd,p(e,:));
            dgradx=(feval(fd,[p(e,1)+deps,p(e,2)])-d)/deps;
            dgrady=(feval(fd,[p(e,1),p(e,2)+deps])-d)/deps;
            dgrad2=dgradx.^2+dgrady.^2;
            dgrad2(dgrad2==0)=1;
            p(e,:)=p(e,:)-[d.*dgradx./dgrad2,d.*dgrady./dgrad2];
        end
        dgnodes(:,:,it) = p;
      end
  else
      for it=tc'
        p = dgnodes(:,:,it);
        deps=sqrt(eps)*max(max(p)-min(p));
        ed = find(mesh.f(abs(mesh.t2f(it,:)),end)<0);
        for id=ed'        
            e = find(philocal(:,id) < 1.e-6);                  
            %d=feval(fd,p(e,:),varargin{:});        
            %dgradx=(feval(fd,[p(e,1)+deps,p(e,2)],varargin{:})-d)/deps;
            %dgrady=(feval(fd,[p(e,1),p(e,2)+deps],varargin{:})-d)/deps;
            d=feval(fd,p(e,:));
            dgradx=(feval(fd,[p(e,1)+deps,p(e,2),p(e,3)])-d)/deps;
            dgrady=(feval(fd,[p(e,1),p(e,2)+deps,p(e,3)])-d)/deps;
            dgradz=(feval(fd,[p(e,1),p(e,2),p(e,3)+deps])-d)/deps;
            dgrad2=dgradx.^2+dgrady.^2+dgradz.^2;
            dgrad2(dgrad2==0)=1;
            p(e,:)=p(e,:)-[d.*dgradx./dgrad2,d.*dgrady./dgrad2,d.*dgradz./dgrad2];
        end
        dgnodes(:,:,it) = p;
      end
  end
end
% if ~isempty(fd) && mesh.porder>1
%   tc=find(mesh.tcurved);
%   for it=tc'
%     p = dgnodes(:,:,it);
%     deps=sqrt(eps)*max(max(p)-min(p));
%     ed = find(mesh.f(abs(mesh.t2f(it,:)),4)<0);
%     for id=ed'        
%         e = find(philocal(:,id) < 1.e-6);        
%         %d=feval(fd,p(e,:),varargin{:});        
%         %dgradx=(feval(fd,[p(e,1)+deps,p(e,2)],varargin{:})-d)/deps;
%         %dgrady=(feval(fd,[p(e,1),p(e,2)+deps],varargin{:})-d)/deps;
%         d=feval(fd,p(e,:));
%         dgradx=(feval(fd,[p(e,1)+deps,p(e,2)])-d)/deps;
%         dgrady=(feval(fd,[p(e,1),p(e,2)+deps])-d)/deps;
%         dgrad2=dgradx.^2+dgrady.^2;
%         dgrad2(dgrad2==0)=1;
%         p(e,:)=p(e,:)-[d.*dgradx./dgrad2,d.*dgrady./dgrad2];
%     end
%     dgnodes(:,:,it) = p;
%   end
% end