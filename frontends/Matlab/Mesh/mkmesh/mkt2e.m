function t2e = mkt2e(t,elemtype,level)
% find neigboring elements depending on level

[nt,ne] = size(t);

t2t = mkt2t(t,elemtype);
list = zeros(1,100);

t2e = zeros(nt,ne);
for j = 1:nt    
  e = t2t(j,:);
  e = e(e>0);  
  k = length(e);
  list(1:k) = e;
  for i = 2:level          
    a = t2t(e,:);
    b = setdiff(a,[j list(1:k)]);
    b = b(b>0);
    m = length(b);
    list((k+1):(k+m)) = b;
    k = k + m;
  end
  t2e(j,1:(k+2)) = [(k+1) j list(1:k)];
end

