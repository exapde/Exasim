function [p1,p2] = fixp(p1,p2,opts,tol,expression)

if nargin<4
    tol = 1e-6;
end
tol2 = tol*tol;

if nargin<5
  expression = @(p) true;
end

dgnodes = 0;
if (size(p1,3) > 1) && (size(p2,3) > 1)
  dgnodes = 1;
  [npe1, nd, ne1] = size(p1);
  p1 = reshape(permute(p1, [2 1 3]), [nd, npe1*ne1]);
  [npe2, nd, ne2] = size(p2);
  p2 = reshape(permute(p2, [2 1 3]), [nd, npe2*ne2]);
end

[nd, np2] = size(p2);

for i = 1:np2  
  if expression(p2(:,i)) 
    e = (p1(1,:) - p2(1,i)).^2;
    for d=2:nd 
        e = e + (p1(d,:)- p2(d,i)).^2;
    end
    [tm,k] = min(e);
    if tm<tol2       
      if (opts==0)
        pa = 0.5*(p2(:,i) + p1(:,k));
        p1(:,k) = pa;
        p2(:,i) = pa;      
      elseif (opts==1)
        p1(:,k) = p2(:,i);            
      elseif (opts==2)
        p2(:,i) = p1(:,k);               
      end
    end
  end
end

if dgnodes==1
  p1 = permute(reshape(p1, [nd npe1 ne1]), [2 1 3]);
  p2 = permute(reshape(p2, [nd npe2 ne2]), [2 1 3]);
end


