function UDG = putsolutiononboundary(UDG, ub, f, perm, ib)

[indf, inde] = find(f == ib);

nfb = length(inde);
ncu = size(ub,2);
for i = 1:nfb
  UDG(perm(:,indf(i)),1:ncu,inde(i)) = ub(:,:,i);  
end

