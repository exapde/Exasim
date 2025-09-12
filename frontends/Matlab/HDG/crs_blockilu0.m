function val = crs_blockilu0(row_ptr, col_ind, val)
% BLOCK_ILU0   In-place block ILU(0) for a CRS matrix with K×K blocks
%
%   val = BLOCK_ILU0(row_ptr, col_ind, val) performs an incomplete LU
%   factorization with zero fill-in (ILU(0)) on a block?sparse matrix
%   stored in CRS.  The input ?val? is modified in place to contain the
%   L (unit?diagonal) and U factors.
%
%   Inputs:
%     row_ptr : (Nnodes+1)x1 array of row?pointers into col_ind/val
%     col_ind : nnz×1 array of column indices (node indices) for each block
%     val     : K×K×nnz array of block entries
%   Output:
%     val     : overwritten with L and U such that L*U is original A

nb = size(val,3);
Nnodes = numel(row_ptr)-1;
%K      = size(val,1);

for i = 1:Nnodes
  for n = 1:nb
    %--- locate and invert U_{ii}
    r1 = row_ptr(i)+1;
    r2 = row_ptr(i+1);
    diag_idx = 0;
    for p = r1:r2
      if col_ind(p)==i
        diag_idx = p;
        break
      end
    end
    if diag_idx==0
      error('ILU0: missing diagonal block at row %d', i);
    end
    Uii_inv = inv(val(:,:,n,diag_idx));
    val(:,:,n,diag_idx) = Uii_inv;        % store U_{ii}^{-1}

    %--- for each neighbor j>i in row i, form L_{ji} and update U_{j,*}
    for p = r1:r2
      j = col_ind(p);
      if j<=i, continue, end

      % find A_{ji} in row j
      rj1 = row_ptr(j)+1;
      rj2 = row_ptr(j+1);
      idx_ji = 0;
      for q = rj1:rj2
        if col_ind(q)==i
          idx_ji = q;
          break
        end
      end
      if idx_ji==0
        continue   % no connection j?i
      end

      % compute L_{ji} = A_{ji} * U_{ii}^{-1}
      Lji = val(:,:,n,idx_ji) * Uii_inv;
      val(:,:,n,idx_ji) = Lji;

      % update U_{j,ell} for each ell>i in the neighbor list of row i
      for pp = r1:r2
        ell = col_ind(pp);
        if ell<=i, continue, end

        % find U_{j,ell} in row j
        idx_jl = 0;
        for qq = rj1:rj2
          if col_ind(qq)==ell
            idx_jl = qq;
            break
          end
        end
        if idx_jl==0
          continue  % ILU(0): skip fill-in
        end

        % Schur complement update: U_{j,ell} -= L_{ji} * U_{i,ell}
        val(:,:,n,idx_jl) = val(:,:,n,idx_jl) - Lji * val(:,:,n,pp);
      end
    end
  end
end

end
