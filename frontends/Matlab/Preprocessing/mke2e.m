function e2e = mke2e(f2e, e2f)

[nfe, ne] = size(e2f);
e2e = zeros(nfe, ne);

for i = 1:ne
  for j = 1:nfe
    k = e2f(j,i); 
    e1 = f2e(1,k);    
    e2 = f2e(3,k);    
    if e1 == i
      e2e(j,i) = e2;      
    elseif e2 == i
      e2e(j,i) = e1;
    else
      error("something wrong");
    end    
  end
end
