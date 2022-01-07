function qdg = graduface(shapft, dgnodes, udg, uh, t2f, perm)

[npf,nfe] = size(perm);
[npe,ncu,ne] = size(udg);
[ngf,npf,nd] = size(shapft);

dgnodes = reshape(dgnodes(perm,:,:),[npf nfe nd ne]);
dgnodes = permute(dgnodes,[1 2 4 3]);

udg = reshape(udg(perm,:,:),[npf nfe ncu ne]);
udg = permute(udg,[1 2 4 3]);
udg = shapft(:,:,1)*reshape(udg,[npf nfe*ne*ncu]);

uh = reshape(uh(:,:,t2f),[npf ncu nfe ne]);
uh = permute(uh,[1 3 4 2]);
uh = shapft(:,:,1)*reshape(uh,[npf nfe*ne*ncu]);

dshapft = reshape(permute(shapft(:,:,2:nd),[1 3 2]),[ngf*(nd-1) npf]);
dpg = dshapft*reshape(dgnodes,[npf nfe*ne*nd]);
dpg = reshape(dpg,[ngf nd-1 nfe ne nd]);
[nlg, jac] = facegeom(dpg,nd);
nlg = reshape(nlg, [ngf*nfe*ne nd]);
jac = reshape(jac, [ngf*nfe*ne 1]);

% ngf*nfe*ncu*ne
% ngf*nfe*nd*ne
% ngf*nfe*ne
% ngf*nfe*ncu*nd*ne
% npf*nfe*ncu*nd*ne
tmp = zeros(ngf*nfe*ne,ncu,nd);
for i = 1:ncu
    for j = 1:nd
        tmp(:,i,j) = (udg(:,i)-uh(:,i)).*nlg(:,j).*jac;
    end
end
tmn = (shapft(:,:,1)')*reshape(tmp,[ngf nfe*ne*ncu*nd]);
tmn = reshape(tmn,[npf nfe ne ncu*nd]);

qdg = zeros(npe, ne, ncu*nd);
for j=1:nfe
    qdg(perm(:,j),:,:) = qdg(perm(:,j),:,:) + reshape(tmn(:,j,:,:),[npf ne ncu*nd]);    
end            


% dpg = permute(reshape(dpg,[ngf nd-1 nfe nd ne]), [1 3 4 2]);    
% dpg = reshape(dpg,[ngf*nfe,nd,nd-1]); 
% [nlg, jacf] = facegeom(dpg,nd);


function [nlg, jac] = facegeom(dpg,nd)

switch nd
    case 1
        jac = 1;
        nlg = [1; -1];        
        nlg = nlg(:);
    case 2
        jac = sqrt(dpg(:,1).^2+dpg(:,2).^2);
        nlg   = [dpg(:,2),-dpg(:,1)];
        nlg   = bsxfun(@rdivide, nlg, jac);
    case 3
        nlg(:,1) = dpg(:,2,1).*dpg(:,3,2) - dpg(:,3,1).*dpg(:,2,2);
        nlg(:,2) = dpg(:,3,1).*dpg(:,1,2) - dpg(:,1,1).*dpg(:,3,2);
        nlg(:,3) = dpg(:,1,1).*dpg(:,2,2) - dpg(:,2,1).*dpg(:,1,2);
        jac = sqrt(nlg(:,1).^2+nlg(:,2).^2+nlg(:,3).^2);
        nlg   = bsxfun(@rdivide, nlg, jac);
    otherwise
        error('Dimension is not implemented');
end







