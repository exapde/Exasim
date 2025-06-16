function uh = faceinsert(uh, uh1, face)

[neb, nfeb] = size(face);
for i = 1:nfeb
  for n = 1:neb
    f = face(n, i);
    uh(:,f) = uh(:,f) + uh1(:,n,i);
  end
end

