function e2f = mke2f(f2e)

ne1 = max(f2e(1,:));
ne2 = max(f2e(3,:));
ne = max(ne1, ne2); 

nfe1 = max(f2e(2,:));
nfe2 = max(f2e(4,:));
nfe = max(nfe1, nfe2); 

e2f = zeros(nfe, ne);
nf = size(f2e, 2);
for i = 1:nf
  e1 = f2e(1,i);
  l1 = f2e(2,i);
  e2 = f2e(3,i);
  e2f(l1, e1) = i;  
  if e2>0
    l2 = f2e(4,i);
    e2f(l2, e2) = i;
  end
end

