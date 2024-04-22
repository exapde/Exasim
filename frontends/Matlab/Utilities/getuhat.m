function uh = getuhat(udg, f2e, perm, ncu)

npf = size(perm,1);
nf = size(f2e,2);
uh = zeros(npf, ncu, nf);
for i = 1:nf
  e1 = f2e(1,i);
  l1 = f2e(2,i);
  uh(:,:,i) = udg(perm(:,l1), 1:ncu, e1);
end
uh = permute(uh, [2 1 3]);
uh = reshape(uh, [ncu npf*nf]);


