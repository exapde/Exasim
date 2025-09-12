function y = pathsolve(A, B, C, D, x)

n = size(A,1);
nep = size(A,3);
m = size(D, 1)/2;

z = x((n*nep+1):end);
for i = 1:nep  
  z(m*(i-1)+1:m*(i+1)) = z(m*(i-1)+1:m*(i+1)) - C(:,:,i)*x((n*(i-1)+1):n*i);
end

ind1 = 1:m;
ind2 = (m+1):2*m;
for i = 1:nep
  z(m*i+1:m*(i+1)) = z(m*i+1:m*(i+1)) - D(ind2,ind1,i)*z(m*(i-1)+1:m*i);  
end

z(m*nep+1:m*(nep+1)) = D(ind2,ind2,nep)*z(m*nep+1:m*(nep+1));
for i = nep:-1:1
  z(m*(i-1)+1:m*i) = D(ind1,ind1,i)*(z(m*(i-1)+1:m*i) - D(ind1,ind2,i)*z(m*i+1:m*(i+1)));
end

y = zeros(n*nep,1);
for i = 1:nep  
  y((n*(i-1)+1):n*i) = A(:,:,i)*(x((n*(i-1)+1):n*i) - B(:,:,i)*z(m*(i-1)+1:m*(i+1)));
end
y = [y; z];

