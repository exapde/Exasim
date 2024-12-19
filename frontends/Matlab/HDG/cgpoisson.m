function [u,K,M,F] = cgpoisson(mesh,master,ubar,par)
mesh.p2 = mesh.p2';
mesh.t2 = mesh.t2';

nn = size(mesh.p2,1);
nd = size(mesh.p2,2);
ne = size(mesh.t2,1);
nv = size(mesh.t2,2);
ng = master.ngv;

il = zeros(nv,nv,ne);
jl = zeros(nv,nv,ne);
for i=1:ne    
    con = mesh.t2(i,:)';    
    com = repmat(con,[1 nv]);%[con con con con con con];
    il(:,:,i) = com;
    jl(:,:,i) = com';        
end

p = mesh.p2(mesh.t2',:);
p = reshape(p,nv,ne,nd);
p = permute(p,[1 3 2]);

if ndims(ubar)<3
    ubar = ubar(mesh.t2');
    ubar = reshape(ubar,nv,ne,1);
    ubar = permute(ubar,[1 3 2]);
end

shapvt  = master.shapvl(:,:,1)'; % ng * nv
dshapvtxi = master.shapvl(:,:,2)';
dshapvtet = master.shapvl(:,:,3)';
dshapvt = reshape(permute(master.shapvl(:,:,2:nd+1),[2 3 1]),[ng*nd nv]);
Ke = zeros(nv,nv,ne);
Me = zeros(nv,nv,ne);
Fe = zeros(nv,ne);
for i=1:ne    
    xg = shapvt*p(:,:,i);
    Jg = reshape(dshapvt*p(:,:,i),[ng,nd,nd]);
    jac = Jg(:,1,1).*Jg(:,2,2) - Jg(:,1,2).*Jg(:,2,1);        
    dshapdx = bsxfun(@times,dshapvtxi,Jg(:,2,2)./jac)-bsxfun(@times,dshapvtet,Jg(:,1,2)./jac);
    dshapdy = bsxfun(@times,dshapvtet,Jg(:,1,1)./jac)-bsxfun(@times,dshapvtxi,Jg(:,2,1)./jac);    
    Ke(:,:,i) = dshapdx'*diag(master.gwvl.*jac)*(dshapdx)+...
                dshapdy'*diag(master.gwvl.*jac)*(dshapdy);    
    Me(:,:,i) = shapvt'*diag(master.gwvl.*jac)*(shapvt);                        
    
    % x = xg(:,1); y = xg(:,2);            
    fg = shapvt*ubar(:,1,i);        
    Fe(:,i) = shapvt'*(master.gwvl.*fg.*jac);        
end

K = sparse(reshape(il,nv*nv*ne,1),reshape(jl,nv*nv*ne,1),reshape(Ke,nv*nv*ne,1));        
M = sparse(reshape(il,nv*nv*ne,1),reshape(jl,nv*nv*ne,1),reshape(Me,nv*nv*ne,1));        
F = sparse(reshape(il(:,1,:),nv*ne,1),ones(nv*ne,1),reshape(Fe,nv*ne,1));    

% u = zeros(nn,1);
% A = par(1)*K+par(2)*M;
% u = A\F;

K(mesh.ib,:) = []; K(:,mesh.ib) = [];
M(mesh.ib,:) = []; M(:,mesh.ib) = [];
F(mesh.ib,:) = [];

u = zeros(nn,1);
A = par(1)*K+par(2)*M;
u(mesh.in) = A\F;

end

