function [L,U] = lu_inplace(A)

[n, m] = size(A);
if n ~= m
    error('Matrix A must be square');
end

U = A;
L = eye(n);

for j = 1:n-1
    for i = j+1:n
        L(i,j) = U(i,j) / U(j,j);
        for k = j+1:n
            U(i,k) = U(i,k) - L(i,j) * U(j,k);
        end
        U(i,j) = 0;  % Optional: zero out for clarity
    end
end

% % LU factorization without pivoting
% % A is an n x n nonsingular matrix
% % LU will contain both L and U:
% %   - U is in the upper triangle (including diagonal)
% %   - L is in the lower triangle (1's are implied on diagonal)
% 
% [n, m] = size(A);
% if n ~= m
%     error('Matrix A must be square.');
% end
% 
% LU = A;
% 
% for j = 1:n-1
%     if LU(j,j) == 0
%         error('Zero pivot encountered at row %d.', j);
%     end
%     for i = j+1:n
%         LU(i,j) = LU(i,j) / LU(j,j);  % Multiplier L(i,j)
%         for k = j+1:n
%             LU(i,k) = LU(i,k) - LU(i,j) * LU(j,k);
%         end
%         LU(i,j+1:n) = LU(i,j+1:n) - LU(i,j) * LU(j,j+1:n);  % Vectorized
%     end
% end
% end
