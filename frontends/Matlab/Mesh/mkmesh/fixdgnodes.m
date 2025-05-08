function dgnodes = fixdgnodes(p,dgnodes,tol)

if nargin<3
    tol = 1e-6;
end

[npe, nd, ne] = size(dgnodes);
if size(p,2)==nd
  p = p';
end

dgnodes = reshape(permute(dgnodes, [2 1 3]), [nd, npe*ne]);

tol2 = tol*tol;
for i = 1:(npe*ne)  
  e = (p(1,:) - dgnodes(1,i)).^2;
  for d=2:nd 
      e = e + (p(d,:)- dgnodes(d,i)).^2;
  end
  [tm,k] = min(e);
  if tm<tol2       
    dgnodes(:,i) = p(:,k);
  end
end

dgnodes = permute(reshape(dgnodes, [nd npe ne]), [2 1 3]);


