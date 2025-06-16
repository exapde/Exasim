function elemblocks = gridpartition3d(Nx, Ny, Nz, Mx, My, Mz, opt)
% Partition a 3D grid into blocks of size Mx x My x Mz

if rem(Nx, Mx) ~= 0
  error("Nx must be divisible by Mx");
end
if rem(Ny, My) ~= 0
  error("Ny must be divisible by My");
end
if rem(Nz, Mz) ~= 0
  error("Nz must be divisible by Mz");
end

nblocks_x = Nx / Mx;
nblocks_y = Ny / My;
nblocks_z = Nz / Mz;

grid = zeros(Nx, Ny, Nz);

if opt == 0
  blocks = zeros(Mx, My, Mz, nblocks_x, nblocks_y, nblocks_z);
  for k = 1:Nz
    for j = 1:Ny
      for i = 1:Nx
        grid(i, j, k) = i + Nx * (j-1) + Nx*Ny*(k-1);
      end
    end
  end

  for kk = 1:nblocks_z
    for jj = 1:nblocks_y
      for ii = 1:nblocks_x
        x_start = (ii-1)*Mx + 1;
        x_end = ii*Mx;
        y_start = (jj-1)*My + 1;
        y_end = jj*My;
        z_start = (kk-1)*Mz + 1;
        z_end = kk*Mz;
        blocks(:,:,:,ii,jj,kk) = grid(x_start:x_end, y_start:y_end, z_start:z_end);
      end
    end
  end

  blocks = reshape(blocks, Mx*My*Mz, nblocks_x*nblocks_y*nblocks_z);

else
  blocks = zeros(Mx, My, Mz, nblocks_z, nblocks_y, nblocks_x);
  for i = 1:Nx
    for j = 1:Ny
      for k = 1:Nz
        grid(i, j, k) = k + Nz*(j-1) + Nz*Ny*(i-1);
      end
    end
  end

  for ii = 1:nblocks_x
    for jj = 1:nblocks_y
      for kk = 1:nblocks_z
        x_start = (ii-1)*Mx + 1;
        x_end = ii*Mx;
        y_start = (jj-1)*My + 1;
        y_end = jj*My;
        z_start = (kk-1)*Mz + 1;
        z_end = kk*Mz;
        blocks(:,:,:,kk,jj,ii) = grid(x_start:x_end, y_start:y_end, z_start:z_end);
      end
    end
  end

  blocks = reshape(blocks, Mx*My*Mz, nblocks_x*nblocks_y*nblocks_z);
end

elemblocks = blocks';

end
