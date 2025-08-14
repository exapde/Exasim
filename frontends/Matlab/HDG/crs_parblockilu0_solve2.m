function x = crs_parblockilu0_solve2(Lind_ji, Uind_ji, Lnum_ji, Unum_ji, val, b)
% BLOCK_ILU0_SOLVE  Solve L*U*x = b using block ILU(0) factors in CRS
%   x = BLOCK_ILU0_SOLVE(row_ptr, col_ind, val, b) performs forward and
%   backward substitution on the block-sparse factors L and U stored in val.

Nnodes = size(Lnum_ji,2);
K = size(val,1);
nb = size(val,3);

% Preallocate y and x
y = zeros(K,nb,Nnodes);
x = zeros(K,nb,Nnodes);

% Forward solve L*y = b (L has unit diagonal)
for i = 1:Nnodes
  yi = b(:,:,i);
  k = Lnum_ji(1,i);
    
  if Lnum_ji(2,i)==1 && k > 0
    p1 = Lind_ji(1,1,i);
    p2 = Lind_ji(k,1,i);
    A = reshape(val(:,:,:,p1:p2),[K K nb*k]);
    B = zeros(K, nb, k);
    for l = 1:k
      j = Lind_ji(l,2,i);
      B(:,:,l) = y(:,:,j);
    end
    B = reshape(B, [K nb*k]);
    C = zeros(K, nb*k);
    for m = 1:nb*k
      C(:,m) = A(:,:,m)*B(:,m);
    end
    C = reshape(C, [K*nb k]);
    yi = yi - reshape(sum(C,2), [K nb]);            
  else
    for l = 1:k
      ptr = Lind_ji(l,1,i);
      j = Lind_ji(l,2,i);
      for n = 1:nb                        
          yi(:,n) = yi(:,n) - val(:,:,n,ptr) * y(:,n,j);
      end
    end
  end
  y(:,:,i) = yi;  
end

% Backward solve U*x = y
for i = Nnodes:-1:1
  yi = y(:,:,i);  
  k = Unum_ji(1,i);  
  rstart = Unum_ji(2,i);  
  if Unum_ji(3,i)==1 && k>0
    p1 = Uind_ji(1,1,i);
    p2 = Uind_ji(k,1,i);
    A = reshape(val(:,:,:,p1:p2),[K K nb*k]);
    B = zeros(K, nb, k);
    for l = 1:k
      j = Uind_ji(l,2,i);
      B(:,:,l) = x(:,:,j);
    end
    B = reshape(B, [K 1 nb*k]);
    C = zeros(K, nb*k);
    for m = 1:nb*k
      C(:,m) = A(:,:,m)*B(:,m);
    end
    C = reshape(C, [K*nb k]);
    yi = yi - reshape(sum(C,2), [K nb]);       
  else    
    for l = 1:k
      ptr = Uind_ji(l,1,i);
      j = Uind_ji(l,2,i);
      for n = 1:nb                        
          yi(:,n) = yi(:,n) - val(:,:,n,ptr) * x(:,n,j);
      end    
    end  
  end
  for n = 1:nb       
    x(:,n,i) = val(:,:,n,rstart) * yi(:,n);
  end    
end

end
