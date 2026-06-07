function u3d = rotate2dfield(u2d, ns)

[npe2d, nc, ne2d] = size(u2d);
npe1d = sqrt(npe2d);

u3d = zeros(npe2d,npe1d,nc,ne2d,ns);

u2d = reshape(u2d, [npe2d 1 nc ne2d]);
for i = 1:ns
  for j = 1:npe1d
    u3d(:,j,:,:,i) = u2d;
  end
end
u3d = reshape(u3d, [npe2d*npe1d nc ne2d*ns]);
