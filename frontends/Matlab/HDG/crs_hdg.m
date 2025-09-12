function [row_ptr, col_ind] = crs_hdg(f2e)

nf = size(f2e,2);

e2f = mke2f(f2e);
[f2f, ~] = mkf2f(f2e, e2f);

row_ptr = zeros(1, nf+1);
for i = 1:nf
  n = 1 + sum(f2f(:,i)>0);
  row_ptr(i+1) = row_ptr(i) + n;
end
col_ind = zeros(1, row_ptr(end));
for i = 1:nf  
  ind = f2f(:,i)>0;
  col_ind((row_ptr(i)+1):row_ptr(i+1)) = [i sort(f2f(ind,i))'];
end

