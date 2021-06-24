function newmesh = refinemesh1(oldmesh,porder,nref)
% refine oldmesh by subdividing each element into nref intervals along each direction 
% and using polymomials of degree porder

check = 0;
ne = oldmesh.ne;
nd = oldmesh.nd;
elemtype = oldmesh.elemtype;
nodetype = oldmesh.nodetype;

% refinement on the master element
[pref,tref] = masternodes(nref,nd,elemtype,0);
npref = size(pref,1);
[ntref,npv]=size(tref);

% shape functions and derivatives on the master element at pref  
shappref = mkshape(oldmesh.porder,oldmesh.plocal,pref,elemtype);
shappref = shappref(:,:,1)';

% new (p,t)
p = zeros(npref, nd, ne);
t = zeros(size(tref,1), size(tref,2), ne);
for i = 1:ne
    p(:,:,i) = shappref*oldmesh.dgnodes(:,1:nd,i);
    t(:,:,i) = tref+(i-1)*npref;
end
if check==1
    p1 = p; 
end
p = reshape(permute(p,[1 3 2]),[npref*ne nd]);
t = reshape(permute(t,[1 3 2]),[size(tref,1)*ne size(tref,2)]);

% remove duplicated nodes
snap = 1e-8;
[~,ix,jx]=unique(round(p/snap)*snap,'rows');
p=p(ix,:);
t=jx(t);

if check==1
    for i=1:ne
        for j = 1:size(tref,1)
            k = (i-1)*size(tref,1) + j;
            if norm(p(t(k,:),:)-p1(tref(j,:),:,i)) > snap
                error('something wrong');
            end
        end    
    end
end

% new refined mesh
newmesh = mkmesh(p,t,porder,oldmesh.bndexpr,elemtype,nodetype);

% higher order master nodes
plocal = masternodes(porder,nd,elemtype,nodetype);
npl=size(plocal,1);

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
    
% dg nodes on the master element
dgref=zeros(npl,nd,ntref);
for dim=1:nd
  for node=1:npv
    dp=philocal(:,node)*pref(tref(:,node),dim)';
    dgref(:,dim,:)=dgref(:,dim,:)+permute(dp,[1,3,2]);
  end
end
dgref = reshape(permute(dgref,[1 3 2]),[npl*ntref nd]);

% shape functions and derivatives on the master element at dgref  
shapdgref = mkshape(oldmesh.porder,oldmesh.plocal,dgref,elemtype);
shapdgref = shapdgref(:,:,1)';

% new dgnodes
newmesh.dgnodes = zeros(npl*ntref, nd, ne);
for i = 1:ne
    newmesh.dgnodes(:,:,i) = shapdgref*oldmesh.dgnodes(:,1:nd,i);
end
newmesh.dgnodes = permute(reshape(newmesh.dgnodes,[npl ntref nd ne]),[1 3 2 4]);
newmesh.dgnodes = reshape(newmesh.dgnodes,[npl nd ntref*ne]);

if check==1
    figure(1); clf; simpplot(pref,tref);
    hold on;
    x = dgref(:,1);
    y = dgref(:,2);
    plot(x(:),y(:),'o');
    
    figure(2); clf; simpplot(p,t); axis on;
    hold on;
    x = newmesh.dgnodes(:,1,:);
    y = newmesh.dgnodes(:,2,:);
    plot(x(:),y(:),'o'); 
end

