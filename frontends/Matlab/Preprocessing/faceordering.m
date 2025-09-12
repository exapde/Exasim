function [ind, nfl, l] = faceordering(f2f)

[p,nf] = size(f2f);
ind = zeros(nf,1);
nfl = zeros(1000,1);

% boundary faces
ind0 = find(f2f(end,:)==0);
ind0 = ind0(:);
nfl(2) = length(ind0);
m = nfl(2);
ind(1:m) = ind0(:);

l = 1;
while (1)
  i1 = nfl(l)+1;
  i2 = nfl(l+1);    
%   ind0 = f2f(:,ind(i1:i2));  % faces connected to the level-l faces
%   ind0 = ind0(ind0>0);       % remove zeros

  % faces connected to the previous level-l faces
  k = 0;
  for i = i1:i2
    for j = 1:p
      tm = f2f(j,ind(i));
      if tm>0
        k = k+1;
        ind0(k) = tm;
      end
    end
  end
  
  if (l==1)
    ind2 = setdiff(ind0(1:k), ind(i1:i2)); % level-(l+1) faces 
  else
    j1 = nfl(l-1)+1;
    j2 = nfl(l);    
    ind1 = union(ind(j1:j2), ind(i1:i2));  % both level-(l-1) and level-l faces   
    ind2 = setdiff(ind0(1:k), ind1);       % level-(l+1) faces     
    if isempty(ind2)
      break;
    end
  end  
  
%   % classify level-(l+1) faces  
%   idx = classifyfaces(ind2, f2f); 
%     
%   % if boundary faces are eliminated, it does not change the connectivities of the faces in ind2.
%   bnd = ind2(idx==1); % boundary faces in ind2
%   cnd = ind2(idx~=1); % interior faces in ind2
%     
%   % reorder and group interior faces to ensure tridiagonal connectivities
%   nbf = findneighborfaces(cnd, f2f);
% %   [cnd(:)'; nbf]
% %   pause
%   
%   l2 = find(nbf(1,:)==2);
%   l3 = nbf(2:3,l2); l3 = unique(l3(:))';
%   for i = 1:length(l3)
%     l3(i) = find(cnd==l3(i));
%   end  
%   l4 = nbf(2:end,l3); l4 = l4(l4 > 0); l4=unique(l4(:)');     
%   l4 = setdiff(l4,union(union(cnd(l2),cnd(l3)),bnd));
%   for i = 1:length(l4)
%     l4(i) = find(cnd==l4(i));
%   end
% 
%   l2 = find(nbf(1,:)==2);
%   l3 = find(nbf(1,:)==3);
%   l4 = find(nbf(1,:)==4);
%   
%   ll = [l2 l3 l4];
%   cnd = cnd(ll);
%   nbf = nbf(:,ll);    
%   ind2 = [bnd; cnd];      % boundary faces and interior faces  
%   
%   ind2 = [9    11    19    27    30    38    53    55   1 7 6 8 5 13 15 16 52 51 50 41 60 57 59 49  3    28    37    56     2     4    17    26    39    48    54  58]';
%   ind2'
%   mbf = findneighborfaces(ind2, f2f);
%   for i = 2:size(mbf,1)
%     for j = 1:size(mbf,2)  
%       if mbf(i,j) > 0
%         m = find(ind2==mbf(i,j));
%         mbf(i,j) = m;
%       end
%     end
%   end  
%   A = eye(length(ind2));
%   for i = 1:length(ind2)
%     for j = 1:mbf(1,i)
%       A(i, mbf(j+1,i)) = rand;
%     end
%   end
%   save A.mat A
%   pause

  i3 = nfl(l+1)+1;
  i4 = nfl(l+1)+length(ind2);     
  ind(i3:i4) = ind2;  
  nfl(l+2) = i4;
  l = l + 1;  
end
nfl = nfl(1:l+1);

end

% function find_many_faces
% 
% end
