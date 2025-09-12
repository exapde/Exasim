function F = assemblelRHS(FE, f2e, elcon)

[npf, nfe, ne] = size(elcon);
ndf = npf*nfe;
nf = size(f2e,2);
ncu = length(FE(:))/(ndf*ne);

FE = reshape(FE, [ncu npf nfe ne]);

F = zeros(ncu, npf, nf);
for i = 1:nf % all faces
  e1 = f2e(1,i);
  l1 = f2e(2,i);
  F(:,:,i) = FE(:,:,l1,e1);
end
for i = 1:nf % interior faces only
  e2 = f2e(3,i);      
  if (e2 > 0)        
    l2 = f2e(4,i);
    ind = elcon(:,l2,e2) - npf*(i-1);    
    F(:,:,i) = F(:,:,i) + FE(:,ind,l2,e2);
  end
end 
F = reshape(F, [ncu, npf*nf]);

