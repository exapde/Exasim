function uhpath = pathextract(uh, fpath, fintf)

nelem = size(fpath,2);
ninte = size(fintf, 1);
nintf = nelem*ninte;
nface = nelem + 1 + nintf;
ncf = size(uh, 1);

uhpath = zeros(ncf, nface);

for i = 1:nelem
  for j = 1:ninte
    f = fintf(j, i);
    uhpath(:,j + ninte*(i-1)) = uh(:,f);
  end
end

for i = 1:nelem
  f = fpath(1,i);
  uhpath(:,nintf + i) = uh(:,f);
end
uhpath(:,end) = uh(:,fpath(2,nelem));


