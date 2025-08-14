function val = crs_parblockilu0(ind_ii, ind_ji, ind_jl, ind_il, num_ji, num_jl, val)
% In-place block ILU(0) for a CRS matrix with K×K blocks

nb = size(val,3);
Nnodes = length(ind_ii);

for i = 1:Nnodes
  diag_idx = ind_ii(i);
  for n = 1:nb    
    val(:,:,n,diag_idx) = inv(val(:,:,n,diag_idx));        
  end
  
  nj = num_ji(i);
  for j = 1:nj
    idx_ji = ind_ji(j,i);
    for n = 1:nb          
      val(:,:,n,idx_ji) = val(:,:,n,idx_ji)*val(:,:,n,diag_idx);        
    end    
    nl = num_jl(j,i);
    for l = 1:nl
      idx_jl = ind_jl(l,j,i);
      idx_il = ind_il(l,j,i);
      for n = 1:nb
        val(:,:,n,idx_jl) = val(:,:,n,idx_jl) - val(:,:,n,idx_ji) * val(:,:,n,idx_il);
      end
    end
  end  
end

end

