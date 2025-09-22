function uh = pathapply(A, B, C, D, b, fpath, fintf, nelems)

ncf = size(b,1);
npaths = length(nelems)-1;

uh = 0*b;
for i = 1:npaths
  ind = (nelems(i)+1):nelems(i+1);  
  x = pathextract(b, fpath(:,ind), fintf(:,ind));  
  y = pathsolve(A(:,:,ind), B(:,:,ind), C(:,:,ind), D(:,:,ind), x(:));     
  uh = pathinsert(uh, reshape(y, ncf, []), fpath(:,ind), fintf(:,ind));
end
