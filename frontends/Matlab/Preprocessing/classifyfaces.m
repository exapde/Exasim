function idx = classifyfaces(ind, f2f)

n = length(ind);
p = size(f2f,1);
f = zeros(p,1);
idx = zeros(n,1);

for i = 1:n              % loop over each face ind(i)  
  j = 0;                 % j is the number of neighbors for face ind(i)
  for k = 1:p            % loop over each neigboring index 
    fk = f2f(k,ind(i));  % get face fk 
    if fk <= 0            
      break;
    elseif ismember(fk, ind) % if face fk is a member of ind
      j = j + 1;
      f(j) = fk;   % f stores neighbors of face ind(i)           
    end
  end  
  
  if (j <= 1)   % face ind(i) has no neighbor or one neighbor
    idx(i) = 1; % set idx(i) = 0 
  else          % face ind(i) has at least two neighbors
    connected_neighbors = 1; % determine if these neighbors are connected 
    for k = 1:j
      for m = 1:j
        if m~=k
          if ismember(f(m),f2f(:,f(k)))==0 % if face f(m) is not a neighbor of face f(k)
            connected_neighbors = 0; % the neighbors are not connected 
            break;
          end
        end
      end
      if (connected_neighbors==0)
        break;
      end
    end
    if connected_neighbors==1 
      idx(i) = 1; % set idx(i) = 1 if the neighbors are connected 
    else
      idx(i) = j; % set idx(i) = j if the neighbors are not connected       
    end
  end
end

end




