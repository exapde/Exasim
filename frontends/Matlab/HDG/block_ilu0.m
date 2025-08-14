function val = block_ilu0(row_ptr, col_ind, val)
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
%     val     : overwritten with L and U such that L*U ? original A

Nnodes = numel(row_ptr)-1;
%K      = size(val,1);

for i = 1:Nnodes
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
  Uii_inv = inv(val(:,:,diag_idx));
  val(:,:,diag_idx) = Uii_inv;        % store U_{ii}^{-1}

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
    Lji = val(:,:,idx_ji) * Uii_inv;
    val(:,:,idx_ji) = Lji;

    % update U_{j,?} for each ?>i in the neighbor list of row i
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

      % Schur?complement update: U_{j,?} -= L_{ji} * U_{i,?}
      val(:,:,idx_jl) = val(:,:,idx_jl) - Lji * val(:,:,pp);
    end
  end
end
end

% 
% function val = block_ilu0(row_ptr, col_ind, val)
% % BLOCK_ILU0 In-place block ILU(0) factorization for block CRS matrices
% %   Performs an incomplete LU factorization with zero fill-in (ILU(0)) on a
% %   matrix stored in compressed row storage (CRS) format with KxK blocks.
% %
% %   INPUTS:
% %     row_ptr : (Nnodes+1)x1 array, row pointers into col_ind and val
% %     col_ind : nnz x1 array, column indices for each block entry
% %     val     : K x K x nnz array of block values (global matrix)
% %   OUTPUT:
% %     val     : overwritten in-place with L and U factors
% %
% %   We compute L_{ji} = A_{ji} * inv(U_{ii}) for j > i, and update
% %   U_{j,ell} = U_{j,ell} - L_{ji} * U_{i,ell} for each ell >= i.
% 
% Nnodes = length(row_ptr) - 1;
% K = size(val,1);
% 
% for i = 1:Nnodes
%     % pointers for row i
%     rstart_i = row_ptr(i) + 1;
%     rend_i   = row_ptr(i+1);
% 
%     % find and invert the diagonal block U_{ii}
%     % locate diagonal block U_{ii} using explicit loop
%     diag_offset = 0;
%     for m = rstart_i:rend_i
%         if col_ind(m) == i
%             diag_offset = m - rstart_i + 1;
%             break;
%         end
%     end
%     if diag_offset == 0
%         error('block_ilu0: missing diagonal block in row %d', i);
%     end
%     
%     diag_idx = rstart_i + diag_offset - 1;
%     Uii_inv = inv(val(:,:,diag_idx));
%     val(:,:,diag_idx) = Uii_inv;
% 
%     % for each neighbor j > i in row i
%     for ptr = rstart_i:rend_i
%         j = col_ind(ptr);
%         if j <= i
%             continue;  % only process lower block rows j>i
%         end
% 
%         % locate block A_{ji} in row j using explicit loop
%         rstart_j = row_ptr(j) + 1;
%         rend_j   = row_ptr(j+1);
%         idx_ji = 0;
%         for m = rstart_j:rend_j
%             if col_ind(m) == i
%                 idx_ji = m;
%                 break;
%             end
%         end
%         if idx_ji == 0
%             continue;  % no connection j <- i
%         end
% 
%         % compute L_{ji} = A_{ji} * U_{ii}^{-1}
%         Lji = val(:,:,idx_ji) * Uii_inv;
%         val(:,:,idx_ji) = Lji;
% 
%         % update U blocks U_{j,ell} for each neighbor ell > i
%         for ptr2 = rstart_i:rend_i
%             ell = col_ind(ptr2);
%             if ell <= i
%                 continue;
%             end
%             
%             % find U_{j,ell} in row j
%             % locate U_{j,ell} entry using explicit loop
%             rel2 = 0;
%             for m2 = rstart_j:rend_j
%                 if col_ind(m2) == ell
%                     rel2 = m2 - rstart_j + 1;
%                     break;
%                 end
%             end
%             if isempty(rel2)
%                 continue;  % ILU(0) no fill-in
%             end
%             idx_jell = rstart_j + rel2 - 1;
% 
%             % U_{j,ell} = U_{j,ell} - L_{ji} * U_{i,ell}
%             val(:,:,idx_jell) = val(:,:,idx_jell) - Lji * val(:,:,ptr2);
%         end
%     end
% end
% end
% 
