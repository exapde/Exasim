function BE = blockjacobi(AE, f2e, elcon)

[npf, nfe, ne] = size(elcon);
ndf = npf*nfe;
nf = size(f2e,2);
n = size(AE,1);
ncu = n/ndf;
ncf = ncu*npf;

AE = reshape(AE,[ncf nfe ncf nfe ne]);

BE = zeros(ncf, ncf, nf);
for i = 1:nf % all faces
  e1 = f2e(1,i);
  l1 = f2e(2,i);
  BE(:,:,i) = reshape(AE(:,l1,:,l1,e1), [ncf ncf]);
end

AE = reshape(AE,[ncu npf nfe ncu npf nfe ne]);
for i = 1:nf % interior faces only
  e2 = f2e(3,i);      
  if (e2 > 0)        
    l2 = f2e(4,i);
    ind = elcon(:,l2,e2) - npf*(i-1);    
    BE(:,:,i) = BE(:,:,i) + reshape(AE(:,ind,l2,:,ind,l2,e2), [ncf ncf]);
  end
end 

for i = 1:nf % all faces
  BE(:,:,i) = inv(BE(:,:,i));
end

