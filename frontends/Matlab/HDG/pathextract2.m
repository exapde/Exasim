function [uh1, uh2] = pathextract2(uh, fpath, fintf, nep)

ne = size(fpath,2);
ninte = size(fintf, 1);
npaths = ne/nep;

ncf = size(uh, 1);
uh1 = zeros(ncf, ninte, npaths, nep);
uh2 = zeros(ncf, npaths, nep+1);

for i = 1:nep
  for n = 1:npaths
    for j = 1:ninte
      f = fintf(j, n + npaths*(i-1));
      uh1(:,j,n,i) = uh(:,f);
    end
    f = fpath(1, n + npaths*(i-1));
    uh2(:,n,i) = uh(:,f);
  end
end
for n = 1:npaths
  f = fpath(2, n + npaths*(nep-1));
  uh2(:,n,nep+1) = uh(:,f);
end

uh1 = reshape(uh1, [ncf*ninte npaths nep]);


