function ub = getuhonboundary(UH, elemcon, f, perm, ib)

[indf, inde] = find(f == ib);

npf = size(perm,1);
nfb = length(inde);
ncu = size(UH,1);

[nfe, ne] = size(f);
elemcon = reshape(elemcon, [npf nfe ne]);

ub = zeros(ncu, npf, nfb);
for i = 1:nfb
  ub(:,:,i) = UH(:,elemcon(:,indf(i),inde(i)));  
end
ub = permute(ub, [2 1 3]);

