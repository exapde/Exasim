function dgnodes = extrudedgnodes(xdg2d,dz)

np2d = size(xdg2d,1);
ne2d = size(xdg2d,3);
np1d = size(dz,1);
nz = size(dz,2);

dgnodes = zeros(np2d*np1d,3,ne2d,nz);
for i = 1:nz
    pm = repmat(dz(:,i)',[np2d 1]);
    for j = 1:ne2d
        dgnodes(:,:,j,i) = [repmat(xdg2d(:,1:2,j),[np1d 1]) pm(:)];
    end
end
dgnodes = reshape(dgnodes,[np2d*np1d,3,ne2d*nz]);

