function [pb, indf, inde] = getponboundary(p, t, f, ib)

[indf, inde] = find(f == ib);
nfb = length(inde);

dim = size(p,1);
n = size(t,1);
elemtype = 0;
if (n - 2^dim) == 0
  elemtype = 1;
end

face = getelemface(dim, elemtype);
nvf = size(face,1);

pb = zeros(dim, nvf, nfb);
for i = 1:nfb
  pi = p(:,t(:,inde(i)));
  pb(:,:,i) = pi(:,face(:,indf(i)));    
end
pb = permute(pb, [2 1 3]);

