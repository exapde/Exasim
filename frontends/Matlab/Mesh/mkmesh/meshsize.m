function [he,hf] = meshsize(mesh)

elemtype = mesh.elemtype;
nd = mesh.nd;
nf = mesh.nf;
nvf = nd;               % number of corners of a face
if nd == 3 && elemtype==1  % hex element
    nvf = nd+1;
end

pf = reshape(mesh.p(mesh.f(:,1:nvf)',:),[nvf nf nd]);

if nvf == 2
    v1 = 1;
    v2 = 2;
elseif nvf == 3
    v1 = [1 2 3];
    v2 = [2 3 1];    
elseif nvf == 4
    v1 = [1 2 3 4];
    v2 = [2 3 4 1];    
end

nv = length(v1);
hf = zeros(nv,nf);
for j = 1:nv
    i1 = v1(j);
    i2 = v2(j);
    hf(j,:) = (pf(i1,:,1)-pf(i2,:,1)).^2;
    for k = 2:nd
        hf(j,:) = hf(j,:) + (pf(i1,:,k)-pf(i2,:,k)).^2;
    end    
end
hf = sqrt(hf);

[ne,nve] = size(mesh.t);
pe = reshape(mesh.p(mesh.t',:),[nve ne nd]);

if nve == 2
    v1 = 1;
    v2 = 2;
elseif nve == 3
    v1 = [1 2 3];
    v2 = [2 3 1];    
elseif nve == 4 && elemtype == 1  % quad element
    v1 = [1 2 3 4];
    v2 = [2 3 4 1];    
elseif nve == 4 && elemtype == 0  % tet element
    v1 = [1 2 3 4 4 4];
    v2 = [2 3 1 1 2 3];        
elseif nve == 8
    v1 = [1 2 3 4 5 6 7 8 1 2 3 4];
    v2 = [2 3 4 1 6 7 8 5 5 6 7 8];    
end

nv = length(v1);
he = zeros(nv,ne);
for j = 1:nv
    i1 = v1(j);
    i2 = v2(j);
    he(j,:) = (pe(i1,:,1)-pe(i2,:,1)).^2;
    for k = 2:nd
        he(j,:) = he(j,:) + (pe(i1,:,k)-pe(i2,:,k)).^2;
    end    
end
he = sqrt(he);





