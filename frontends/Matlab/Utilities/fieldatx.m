function [udgx] = fieldatx(mesh,master,UDG,x,nref)

if size(mesh.t,1) < size(mesh.t,2)    
    mesh.t = mesh.t';
    mesh.p = mesh.p';
end
mesh.nd = master.nd;
mesh.nodetype = master.nodetype;
mesh.elemtype = master.elemtype;
mesh.plocal = master.xpe;
mesh.tlocal = master.telem;

nc = size(UDG,2);
nd = mesh.nd;
[ne,nv] = size(mesh.t);

% element centers
p = reshape(mesh.p(mesh.t',:),[nv ne nd]);
xm = reshape(mean(p,1),[ne nd]);

% get fields and coordinates
udg = permute(UDG,[1 3 2]);
pdg = permute(mesh.dgnodes(:,1:nd,:),[1 3 2]);

% subdivision  
[pdg,udg] = scalarrefine(mesh,pdg,udg,nref);              
npe = size(pdg,1);
pdg = reshape(pdg,[npe*ne nd]);
udg = reshape(udg,[npe*ne nc]);
indp = (1:npe)';

nx = size(x,1);
udgx = zeros(nx,nc);
for i = 1:nx
    if (rem(i,1000)==0)
        [i nx]
    end
    distance = sqrt((x(i,1)-xm(:,1)).^2 + (x(i,2)-xm(:,2)).^2);    
    [~,inde] = sort(distance);
    inde = inde(1:min(100, length(inde)));
    ind = (inde(:)'-1)*npe+indp;
    xa = pdg(ind,:);
    distance = sqrt((x(i,1)-xa(:,1)).^2 + (x(i,2)-xa(:,2)).^2);    
    [~,imin] = min(distance);
    udgx(i,:) = udg(ind(imin),:);    
end


function [pref,uref] = scalarrefine(mesh,p,u,nref)

[npl, nt, nd] = size(p);
porder=mesh.porder;
plocal=mesh.plocal;
tlocal=mesh.tlocal;

if isempty(nref), nref=ceil(log2(max(porder,1))); end
if mesh.elemtype==0  
    A0=koornwinder(plocal(:,1:nd),porder);
    [plocal,tlocal]=uniref(plocal,tlocal,nref);    
    A=koornwinder(plocal(:,1:nd),porder)/A0;
else
    A0=tensorproduct(plocal(:,1:nd),porder);
    m = porder*(nref+1)+1;     
    if nd==2
        [plocal,tlocal]=squaremesh(m,m,1);
    else
        [plocal,tlocal]=cubemesh(m,m,m,1);
    end
    A=tensorproduct(plocal(:,1:nd),porder)/A0;  
end

npln=size(plocal,1);
% t = kron(ones(nt,1),tlocal)+kron(npln*(0:nt-1)',0*tlocal+1);
% ne = size(t,1);
% np = npln*nt;

nc = size(u,3);
uref = reshape(A*reshape(u,npl,nt*nc),[npln nt nc]);
pref = reshape(A*reshape(p,npl,nt*nd),[npln nt nd]);

% uref = reshape(A*reshape(u,npl,nt*nc),[np nc]);
% uref = reshape(uref(t',:), [size(t,2) ne nc]);
% 
% pref = reshape(A*reshape(p,npl,nt*nd),[np,nd]);
% pref = reshape(pref(t',:),[size(t,2) ne nd]);

%function plotTetrahedron(x, p)






