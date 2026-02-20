L1 = eigs(K);

sz = ncu * npf;
I = (1:sz);
for i = 1:nf  
  K(I,I) = BE(:,:,i)*K(I,I);
  I = I + sz;
end

L2 = eigs(K);
