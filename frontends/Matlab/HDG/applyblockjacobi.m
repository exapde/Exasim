function w = applyblockjacobi(BE, v)

ncf = size(BE,1);
nf = size(BE,3);
v = reshape(v, [ncf, nf]);
w = zeros(ncf, nf);
for i = 1:nf
  w(:,i) = BE(:,:,i)*v(:,i);
end



