function nbf = findneighborfaces(ind, f2f)

n = length(ind);
p = size(f2f,1);
nbf = zeros(p+1,n);

for i = 1:n              % loop over each face ind(i)  
  j = 0;                 % j is the number of neighbors for face ind(i)
  for k = 1:p            % loop over each neigboring index 
    fk = f2f(k,ind(i));  % get face fk 
    if fk <= 0            
      break;
    elseif ismember(fk, ind) % if face fk is a member of ind
      j = j + 1;
      nbf(j+1,i) = fk;       % stores neighbors of face ind(i)           
    end
  end    
  nbf(1,i) = j;              % store the number of neighbors for face ind(i)
end

end




