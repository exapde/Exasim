function [ind_ii, ind_ji, ind_jl, ind_il, num_ji, num_jl, Lind_ji, Uind_ji, Lnum_ji, Unum_ji] = crs_indexingilu0(row_ptr, col_ind, nfe)

N = numel(row_ptr)-1;
M = 2*(nfe - 1);

ind_ii = zeros(1, N);
ind_ji = zeros(M, N);
ind_jl = zeros(M, M, N);
ind_il = zeros(M, M, N);

num_ji = zeros(1, N);
num_jl = zeros(M, N);

for i = 1:N
  r1 = row_ptr(i)+1;
  r2 = row_ptr(i+1);
  diag_idx = 0;
  for p = r1:r2
    if col_ind(p)==i
      diag_idx = p;
      break
    end
  end
  if diag_idx==0
    error('ILU0: missing diagonal block at row %d', i);
  end

  ind_ii(i) = diag_idx;
  
  k = 0;  
  for p = r1:r2
    j = col_ind(p);
    if j<=i, continue, end

    rj1 = row_ptr(j)+1;
    rj2 = row_ptr(j+1);
    idx_ji = 0;
    for q = rj1:rj2
      if col_ind(q)==i
        idx_ji = q;
        break
      end
    end
    if idx_ji==0
      continue   
    end
    
    k = k+1;
    ind_ji(k, i) = idx_ji;
    num_ji(i) = k;
    
    m = 0;
    for pp = r1:r2
      ell = col_ind(pp);
      if ell<=i, continue, end

      idx_jl = 0;
      for qq = rj1:rj2
        if col_ind(qq)==ell
          idx_jl = qq;
          break
        end
      end
      if idx_jl==0
        continue  % ILU(0): skip fill-in
      end
      
      m = m+1;
      ind_jl(m, k, i) = idx_jl;
      ind_il(m, k, i) = pp;
      num_jl(k,i) = m;
    end
  end
end

Lind_ji = zeros(M,2,N);
Lnum_ji = zeros(2,N);
for i = 1:N
  rstart = row_ptr(i) + 1;
  rend   = row_ptr(i+1);  
  k = 0;
  for ptr = rstart:rend
      j = col_ind(ptr);
      if j < i
        k = k+1;
        Lind_ji(k,1,i) = ptr;                        
        Lind_ji(k,2,i) = j;           
      end
  end  
  Lnum_ji(1,i) = k;
  Lnum_ji(2,i) = 1;
  for l = 2:k
    if Lind_ji(l,1,i)-Lind_ji(l-1,1,i) ~=1
      Lnum_ji(2,i) = 0;
      break;
    end
  end  
end

Uind_ji = zeros(M,2,N);
Unum_ji = zeros(3,N);
for i = N:-1:1
  rstart = row_ptr(i) + 1;
  rend   = row_ptr(i+1);
  k = 0;
  for ptr = rstart:rend
      ell = col_ind(ptr);
      if ell > i
        k = k+1;
        Uind_ji(k,1,i) = ptr;                        
        Uind_ji(k,2,i) = ell;                        
      end
  end  
  Unum_ji(1,i) = k;
  Unum_ji(2,i) = rstart;
  Unum_ji(3,i) = 1;
  for l = 2:k
    if Uind_ji(l,1,i)-Uind_ji(l-1,1,i) ~=1
      Unum_ji(i) = 0;
      break;
    end
  end  
end

end
