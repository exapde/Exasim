function uh1 = faceextract(uh, face)

[neb, nfeb] = size(face);

ncf = size(uh, 1);
uh1 = zeros(ncf, neb, nfeb);

for i = 1:nfeb
  for n = 1:neb
    f = face(n, i);
    uh1(:,n,i) = uh(:,f);
  end
end

