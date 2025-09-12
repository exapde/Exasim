function test_block_ilu0()
% TEST_BLOCK_ILU0  Unit test for block_ilu0 with a small random block-sparse matrix

% Parameters
K = 2;           % block size
Nnodes = 4;      % number of nodes (rows/cols)
% Construct a simple banded sparsity pattern with diagonal + first off-diagonals
max_band = 1;
% Preallocate connectivity in CRS form
dummy_val = [];
row_ptr = zeros(Nnodes+1,1);
col_ind = [];
val = [];
row_ptr(1) = 0;

% Build row_ptr and col_ind for banded structure
total_nnz = 0;
for i = 1:Nnodes
    cols = max(1,i-max_band):min(Nnodes,i+max_band);
    nnz_row = numel(cols);
    col_ind = [col_ind; cols'];
    total_nnz = total_nnz + nnz_row;
    row_ptr(i+1) = total_nnz;
end

% Initialize random block values ensuring diagonal blocks are invertible
val = zeros(K,K,total_nnz);
idx = 1;
for i = 1:Nnodes
    cols = col_ind(row_ptr(i)+1:row_ptr(i+1));
    for j = cols'
        block = randn(K);
        if j == i
            % Make diagonal block diagonally dominant for invertibility
            block = block + K*eye(K);
        end
        val(:,:,idx) = block;
        idx = idx + 1;
    end
end

% Keep original full matrix
A_full = zeros(Nnodes*K);
idx = 1;
for i = 1:Nnodes
    for ptr = row_ptr(i)+1:row_ptr(i+1)
        j = col_ind(ptr);
        rows = (i-1)*K + (1:K);
        cols = (j-1)*K + (1:K);
        A_full(rows,cols) = val(:,:,idx);
        idx = idx + 1;
    end
end

% Perform in-place ILU(0)
val_fact = block_ilu0(row_ptr, col_ind, val);

% Extract L and U into full matrices
L_full = eye(Nnodes*K);
U_full = zeros(Nnodes*K);
idx = 1;
for i = 1:Nnodes
    for ptr = row_ptr(i)+1:row_ptr(i+1)
        j = col_ind(ptr);
        rows = (i-1)*K + (1:K);
        cols = (j-1)*K + (1:K);
        block = val_fact(:,:,idx);
        if j < i
            % L_{ij}
            L_full(rows,cols) = block;
        elseif j == i
            % diag block holds inv(U_ii)
            Uii_inv = block;
            Uii = inv(Uii_inv);
            U_full(rows,cols) = Uii;
        else
            % U_{ij}
            U_full(rows,cols) = block;
        end
        idx = idx + 1;
    end
end

% Reconstruct A and compare
A_rec = L_full * U_full;
diff_norm = norm(A_rec - A_full, 'fro');

fprintf('Block ILU(0) Factorization Test: Frobenius norm ||L*U - A|| = %g\n', diff_norm);
if diff_norm < 1e-10
    fprintf('Block ILU(0) Factorization Test: PASSED.\n');
else
    fprintf('Block ILU(0) Factorization Test: FAILED.\n');
end

% --- Now test the ILU0 solver ---
% Generate random right-hand side
b_rand = rand(Nnodes*K,1);
% Solve L U x = b using the block ILU0 factors
x_sol = block_ilu0_solve(row_ptr, col_ind, val_fact, b_rand);
% Compute residual norm ||A*x - b||
res = norm(A_full*x_sol - b_rand, 'fro');
fprintf('Block ILU(0) Solver test: ||A*x - b|| = %g\n', res);
if res < 1e-10
    fprintf('Block ILU(0) Solver Test: PASSED.\n');
else
    fprintf('Block ILU(0) Solver Test: FAILED.\n');
end

end
