function x = crs_parblockilu0_solve(row_ptr, col_ind, val, b)
% BLOCK_ILU0_SOLVE  Solve L*U*x = b using block ILU(0) factors in CRS
%   x = BLOCK_ILU0_SOLVE(row_ptr, col_ind, val, b) performs forward and
%   backward substitution on the block-sparse factors L and U stored in val.
%
%   Inputs:
%     row_ptr : (Nnodes+1)x1 row pointers
%     col_ind : nnz x1 column indices for each block
%     val     : K x K x nnz array, containing inverted U_ii on diagonals,
%               L_{ji} for j>i in lower positions, and U_{i?} for ?>=i in upper.
%     b       : (K*Nnodes)x1 right-hand side vector
%   Output:
%     x       : (K*Nnodes)x1 solution vector

Nnodes = length(row_ptr) - 1;
K = size(val,1);
nb = size(val,3);

% Preallocate y and x
y = zeros(K,nb,Nnodes);
x = zeros(K,nb,Nnodes);

% Forward solve L*y = b (L has unit diagonal)
for i = 1:Nnodes
  rstart = row_ptr(i) + 1;
  rend   = row_ptr(i+1);  
  yi = b(:,:,i);
  for ptr = rstart:rend
      j = col_ind(ptr);
      if j < i
        for n = 1:nb                        
            yi(:,n) = yi(:,n) - val(:,:,n,ptr) * y(:,n,j);
        end
      end
  end
  y(:,:,i) = yi;  
end

% Backward solve U*x = y
for i = Nnodes:-1:1
  rstart = row_ptr(i) + 1;
  rend   = row_ptr(i+1);
  yi = y(:,:,i);
  
  for ptr = rstart:rend
      ell = col_ind(ptr);
      if ell > i
        for n = 1:nb                        
            yi(:,n) = yi(:,n) - val(:,:,n,ptr) * x(:,n,ell);
        end
      end
  end
  
  for n = 1:nb       
    x(:,n,i) = val(:,:,n,rstart) * yi(:,n);
  end    
end

end
