function uh = pathinsert(uh, uhpath, fpath, fintf)

nelem = size(fpath,2);
ninte = size(fintf, 1);
nintf = nelem*ninte;

for i = 1:nelem
  for j = 1:ninte
    f = fintf(j, i);
    uh(:,f) = uh(:,f) + uhpath(:,j + ninte*(i-1));
  end
end

for i = 1:nelem
  f = fpath(1,i);
  uh(:,f) = uh(:,f) + uhpath(:,nintf + i);
end
f = fpath(2,nelem);
uh(:,f) = uh(:,f) + uhpath(:,end);



