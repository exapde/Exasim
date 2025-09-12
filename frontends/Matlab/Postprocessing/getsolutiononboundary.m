function ub = getsolutiononboundary(UDG, f, perm, ib)

[indf, inde] = find(f == ib);

npf = size(perm,1);
nfb = length(inde);
ncu = size(UDG,2);

ub = zeros(npf, ncu, nfb);
for i = 1:nfb
  ub(:,:,i) = UDG(perm(:,indf(i)),:,inde(i));  
end
