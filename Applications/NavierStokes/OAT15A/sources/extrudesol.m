function UDG3D = extrudesol(UDG2D, porder, nz)
[np2d,nc,ne2d] = size(UDG2D);
np1d = porder+1;
UDG3D = zeros(np2d*np1d,nc,ne2d,nz);
for j = 1:ne2d
    UDG2Di = repmat(UDG2D(:,:,j),[np1d 1]);
    for i = 1:nz
        for k = 1:nc
            UDG3D(:,k,j,i) = UDG2Di(:,k); % r
        end
    end
end
UDG3D = reshape(UDG3D,[np2d*np1d,nc,ne2d*nz]);
