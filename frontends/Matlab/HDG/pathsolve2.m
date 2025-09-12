function [x1, x2] = pathsolve2(A, B1, B2, C1, C2, D1, D2, DL, DU, x1, x2)

% A : q * q * npaths * nep
% B1 : q * m * npaths * nep
% B2 : q * m * npaths * nep
% C1 : m * q * npaths * nep
% C2 : m * q * npaths * nep
% D1 : m * m * npaths * nep
% D2 : m * m * npaths * nep
% DL : m * m * npaths * nep
% DU : m * m * npaths * nep

% x1 : q * npaths * nep 
% x2 : m * npaths * (nep+1)

npaths = size(A,3);
nep = size(A,4);

for i = 1:nep  
  for n = 1:npaths
    x2(:,n,i) = x2(:,n,i) - C1(:,:,n,i)*x1(:,n,i);
  end  
  for n = 1:npaths
    x2(:,n,i+1) = x2(:,n,i+1) - C2(:,:,n,i)*x1(:,n,i);
  end
  
  for n = 1:npaths
    x2(:,n,i+1) = x2(:,n,i+1) - DL(:,:,n,i)*x2(:,n,i);  
  end  
end

for n = 1:npaths
  x2(:,n,nep+1) = D2(:,:,n,nep)*x2(:,n,nep+1);
end
for i = nep:-1:1
  for n = 1:npaths
    x2(:,n,i) = D1(:,:,n,i)*(x2(:,n,i) - DU(:,:,n,i)*x2(:,n,i+1));
  end
end

for i = 1:nep  
  for n = 1:npaths
    x1(:,n,i) = A(:,:,n,i)*(x1(:,n,i) - B1(:,:,n,i)*x2(:,n,i) - B2(:,:,n,i)*x2(:,n,i+1));
  end
end

