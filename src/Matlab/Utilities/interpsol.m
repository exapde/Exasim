function [UDG2] = interpsol(mesh2, mesh1, UDG1, nref)

npe = size(mesh2.dgnodes,1);
ne = size(mesh2.dgnodes,3);
ncu = size(UDG1,2);

% get CG field on mesh1
[~,cgelcon,rowent2elem,~,cgent2dgent] = mkcgent2dgent(mesh1.dgnodes(:,1:mesh1.nd,:),1e-5);
UCG1 = dg2cg(UDG1, cgelcon, cgent2dgent, rowent2elem);

for i = 1:ncu
    figure(i); clf; scaplot(mesh1,UCG1(:,i,:),[],1); axis tight; axis off; colormap jet;
end

% get CG field on mesh2
[cgnodes,cgelcon] = mkcgent2dgent(mesh2.dgnodes(:,1:mesh2.nd,:),1e-6);
UCG2 = fieldatx(mesh1,UCG1,cgnodes,nref);

UDG2 = zeros(npe,ncu,ne);
for i = 1:ne
    UDG2(:,:,i) = UCG2(cgelcon(:,i),:);
end

for i = 1:ncu
    figure(i+ncu); clf; scaplot(mesh2,UDG2(:,i,:),[],1); axis tight; axis off; colormap jet;
end




