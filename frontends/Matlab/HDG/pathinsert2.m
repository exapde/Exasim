function uh = pathinsert2(uh, uh1, uh2, fpath, fintf, nep)

ne = size(fpath,2);
ninte = size(fintf, 1);
npaths = ne/nep;

ncf = size(uh, 1);
uh1 = reshape(uh1, [ncf ninte npaths nep]);
uh2 = reshape(uh2, [ncf npaths nep+1]);

for i = 1:nep
  for n = 1:npaths
    for j = 1:ninte
      f = fintf(j, n + npaths*(i-1));
      uh(:,f) = uh(:,f) + uh1(:,j,n,i);
    end
    f = fpath(1, n + npaths*(i-1));
    uh(:,f) = uh(:,f) + uh2(:,n,i);
  end
end
for n = 1:npaths
  f = fpath(2, n + npaths*(nep-1));
  uh(:,f) = uh(:,f) + uh2(:,n,nep+1);
end






