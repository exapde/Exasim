function x = block_ilu0_solve(row_ptr, col_ind, val, b)
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

% Preallocate y and x
y = zeros(K*Nnodes,1);
x = zeros(K*Nnodes,1);

% Forward solve L*y = b (L has unit diagonal)
for i = 1:Nnodes
    % extract b_i
    yi = b((i-1)*K + (1:K));
    % row pointers for i
    rstart = row_ptr(i) + 1;
    rend   = row_ptr(i+1);
    % subtract contributions from L_ij * y_j for j<i
    for ptr = rstart:rend
        j = col_ind(ptr);
        if j < i
            Lij = val(:,:,ptr);
            yj = y((j-1)*K + (1:K));
            yi = yi - Lij * yj;
        end
    end
    y((i-1)*K + (1:K)) = yi;
end

% Backward solve U*x = y
for i = Nnodes:-1:1
    % extract y_i
    yi = y((i-1)*K + (1:K));
    rstart = row_ptr(i) + 1;
    rend   = row_ptr(i+1);
    % subtract U_{i?} * x_? for ?>i
    for ptr = rstart:rend
        ell = col_ind(ptr);
        if ell > i
            Uil = val(:,:,ptr);
            xell = x((ell-1)*K + (1:K));
            yi = yi - Uil * xell;
        end
    end
    % apply diagonal inverse U_{ii}^{-1} stored in val
    % find diagonal index
    diag_idx = 0;
    for ptr = rstart:rend
        if col_ind(ptr) == i
            diag_idx = ptr;
            break;
        end
    end
    if diag_idx == 0
        error('block_ilu0_solve: missing diagonal block at row %d', i);
    end
    Uii_inv = val(:,:,diag_idx);
    xi = Uii_inv * yi;
    x((i-1)*K + (1:K)) = xi;
end

end
