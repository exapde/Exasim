function [dgnodes,elemtype,perm] = createdgnodes(p,t,f,curvedboundary,curvedboundaryexpr,porder)
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

[plocal,~,~,~,perm] = masternodes(porder,nd,elemtype);

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

if ~isempty(curvedboundaryexpr) && porder>1 && max(curvedboundary)>0 
    disp('Project dgnodes onto the curved boundaries...');  
    fd = curvedboundaryexpr;
    nd=size(p,1);
    [nfe,ne] = size(f);  
    if nd == 2      
        for i = 1:ne
            for j = 1:nfe
                if f(j,i)~=0 % boundary element
                    k = abs(f(j,i)); % get boundary index
                    if curvedboundary(k)==1 % if this boundary is curved
                        p = dgnodes(perm(:,j),:,i);
                        deps = sqrt(eps)*max(max(p)-min(p));     
                        d = (fd{k}(p'))';
                        dgradx = ((fd{k}([p(:,1)+deps,p(:,2)]'))'-d)/deps;
                        dgrady = ((fd{k}([p(:,1),p(:,2)+deps]'))'-d)/deps;
                        dgrad2 = dgradx.^2+dgrady.^2;
                        dgrad2(dgrad2==0) = 1;
                        p = p-[d.*dgradx./dgrad2,d.*dgrady./dgrad2];                  
                        dgnodes(perm(:,j),:,i) = p;                        
                    end
                end
            end
        end
    elseif nd==3
        for i = 1:ne
            for j = 1:nfe
                if f(j,i)~=0 % boundary element
                    k = abs(f(j,i)); % get boundary index
                    if curvedboundary(k)==1 % if this boundary is curved                
                        p = dgnodes(perm(:,j),:,i);
                        deps = sqrt(eps)*max(max(p)-min(p));     
                        d = (fd{k}(p'))';
                        dgradx = ((fd{k}([p(:,1)+deps,p(:,2),p(:,3)]'))'-d)/deps;
                        dgrady = ((fd{k}([p(:,1),p(:,2)+deps,p(:,3)]'))'-d)/deps;
                        dgradz = ((fd{k}([p(:,1),p(:,2),p(:,3)+deps]'))'-d)/deps;
                        dgrad2 = dgradx.^2+dgrady.^2+dgradz.^2;
                        dgrad2(dgrad2==0) = 1;
                        p = p-[d.*dgradx./dgrad2,d.*dgrady./dgrad2,d.*dgradz./dgrad2];
                        dgnodes(perm(:,j),:,i) = p;         
                    end
                end
            end
        end        
    end        
end

end
