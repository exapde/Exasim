function idx = nodes2indices(nodes, row_ptr, col_ind)

N = length(nodes);
idx = zeros(N,N);
for n = 1:N
  i = nodes(n);
  row_start = row_ptr(i) + 1;
  row_end = row_ptr(i+1);      
  for m = 1:N
    if m == n
      idx(n,m) = row_start;
    else
      for k = row_start+1:row_end
        if col_ind(k) == nodes(m)
          idx(n,m) = k;
          break;
        end        
      end  
    end
  end  
end

end
