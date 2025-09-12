function ind = findinvalidxi(xi, tol, elemtype)

xdim = ndims(xi);
if xdim == 3
  sz = size(xi);
  xi = reshape(permute(xi,[1 3 2]), [sz(1)*sz(3) sz(2)]);
end

nd = size(xi,2);

if nd==2
  ind = find(xi(:,1)<-tol | xi(:,2)<-tol | xi(:,1)>1+tol | xi(:,2)>1+tol); 
  if elemtype==0
    ind1 = find(xi(:,1)+xi(:,2)>1+tol);
    ind = union(ind,ind1);
  end
end

