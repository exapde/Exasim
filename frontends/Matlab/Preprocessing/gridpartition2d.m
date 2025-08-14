function elemblocks = gridpartition2d(Nx, Ny, Mx, My, opt)

if rem(Nx, Mx) ~= 0
  error("Nx must be divisible by Mx");
end
if rem(Ny, My) ~= 0
  error("Ny must be divisible by My");
end
  
nblocks_x = Nx/Mx;
nblocks_y = Ny/My;
grid = zeros(Nx, Ny);

if opt == 0
  blocks = zeros(Mx, My, nblocks_x, nblocks_y);
  for j = 1:Ny
    for i = 1:Nx
      grid(i, j) = i + Nx*(j-1);
    end
  end
  
  for j = 1:nblocks_y
    for i = 1:nblocks_x
      x_start = (i-1)*Mx + 1;
      x_end = i*Mx;
      y_start = (j-1)*My + 1;
      y_end = j*My;
      blocks(:,:,i,j) = grid(x_start:x_end, y_start:y_end);
    end
  end  
else
  blocks = zeros(Mx, My, nblocks_y, nblocks_x);
  for i = 1:Nx
    for j = 1:Ny
      grid(i, j) = j + Ny*(i-1);
    end
  end
  
  for i = 1:nblocks_x
    for j = 1:nblocks_y    
      x_start = (i-1)*Mx + 1;
      x_end = i*Mx;
      y_start = (j-1)*My + 1;
      y_end = j*My;
      blocks(:,:,j,i) = grid(x_start:x_end, y_start:y_end);
    end
  end  
end

elemblocks = reshape(blocks, Mx*My, nblocks_x*nblocks_y)';

end




