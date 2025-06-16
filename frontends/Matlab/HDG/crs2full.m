function A_full = crs2full(row_ptr, col_ind, val)

K = size(val,1);
Nnodes = length(row_ptr)-1;

%  full matrix
A_full = zeros(Nnodes*K, Nnodes*K);
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
