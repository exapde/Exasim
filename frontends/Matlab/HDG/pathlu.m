function [A, C, D] = pathlu(A, B, C, D)

nep = size(A,3);
m = size(D, 1)/2;

for i = 1:nep
  A(:,:,i) = inv(A(:,:,i));    
  C(:,:,i) = C(:,:,i)*A(:,:,i);    
  CB = C(:,:,i)*B(:,:,i);
  D(:,:,i) = D(:,:,i) - CB;
  if i < nep    
    D(1:m,1:m,i+1) = D(1:m,1:m,i+1) - CB(m+1:2*m,m+1:2*m); 
  end
end

% LU decomposition of tridagonal system
for i = 1:nep
  D(1:m,1:m,i) = inv(D(1:m,1:m,i));
  D(m+1:2*m,1:m,i) = D(m+1:2*m,1:m,i)*D(1:m,1:m,i);   
  if i < nep
    D(1:m,1:m,i+1) = D(1:m,1:m,i+1) - D(m+1:2*m,1:m,i)*D(1:m,m+1:2*m,i);
  else
    D(m+1:2*m,m+1:2*m,i) = D(m+1:2*m,m+1:2*m,i) - D(m+1:2*m,1:m,i)*D(1:m,m+1:2*m,i);
    D(m+1:2*m,m+1:2*m,i) = inv(D(m+1:2*m,m+1:2*m,i));
  end
end


