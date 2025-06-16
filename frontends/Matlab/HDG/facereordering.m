function [row_ptr, col_ind, face, f2e, row_ptr2, col_ind2, face2, color, idr1, idr2, idr3, idx1, idx2, idx3] = facereordering(elem, t, e2f, elemtype, nd)

[row_ptr, col_ind, f2e, f2f] = facereorder(elem(1,:), t, elemtype, nd);

nf = size(f2e,2);
nsub = size(elem, 1);
face = zeros(nsub,nf);
for n = 1:nsub    
  for i = 1:nf
    e = f2e(1,i);
    l = f2e(2,i);
    face(n,i) = e2f(l,elem(n,e));
  end
end

ind2 = find(f2e(end,:)>0);
nf2 = length(ind2);
f2e2 = f2e(:,ind2);
[row_ptr2, col_ind2] = crs_hdg(f2e2);

face2 = zeros(nsub,nf2);
for n = 1:nsub    
  for i = 1:nf2
    e = f2e2(1,i);
    l = f2e2(2,i);
    face2(n,i) = e2f(l,elem(n,e));
  end
end

nb = row_ptr(end);
color = zeros(6,nb);
for i = 1:nf % loop over each local face i
  row_start = row_ptr(i) + 1;
  row_end = row_ptr(i+1);      
  for k = row_start:row_end
    j = col_ind(k);  % neighboring face j    
    color(5,k) = i;
    color(6,k) = j;
  end
end

nfe = size(f2f,1)/2 + 1; 
ind1 = find(f2e(end,:)==0);
nf1 = length(ind1);

count1 = 0;
count2 = 0;
count3 = 0;
idr1 = zeros(nfe-1, nf1);
idr2 = zeros(nfe-2, nf1);
idr3 = zeros(nfe-3, nf1);
idx1 = zeros(nfe-1, nfe-1, nf1);
idx2 = zeros(nfe-2, nfe-2, nf1);
idx3 = zeros(nfe-3, nfe-3, nf1);

i = 1;
while (i <= nf1)  
  % find end of this run
  next = i + 1;
  while next <= nf1 && f2e(1,next) == f2e(1,i)
      next = next + 1;
  end
  cnt = next - i;
  if cnt>3
    error("An element on the interface have more than 3 boundary faces.");
  end
      
  row_start = row_ptr(i) + 1;
  row_end = row_ptr(i+1);
  fi = col_ind((row_start+1):row_end);    
  if cnt == 1 % one boundary face on one element
    count1 = count1 + 1;
    % A1 : 1 x  1
    color(1,row_start) = 1; 
    color(2,row_start) = i;
    color(3,row_start) = i;
    color(4,row_start) = i;               
    for n = 1:length(fi) 
      % B1 : 1 x (nfe - 1)    
      k = row_start + n;
      color(1,k) = 4; 
      color(2,k) = i;
      color(3,k) = 1;
      color(4,k) = n;      
      j = fi(n);
      for l = (row_ptr(j) + 1):row_ptr(j+1)
        if col_ind(l) == i                     
          % C1 : (nfe - 1) x 1    
          color(1,l) = 7; 
          color(2,l) = i;
          color(3,l) = n;
          color(4,l) = 1;                
          break;
        end
      end
    end 
    
    idr1(:,count1) = fi - nf1;
    idx1(:,:,count1) = nodes2indices(fi - nf1, row_ptr2, col_ind2);        
  elseif cnt == 2 % two boundary faces on one element
    % A2 : 2 x 2
    count2 = count2 + 1;
    color(1,row_start) = 2;
    color(2,row_start) = count2;
    color(3,row_start) = 1;
    color(4,row_start) = 1;
    
    color(1,row_start+1) = 2;
    color(2,row_start+1) = count2;
    color(3,row_start+1) = 1;
    color(4,row_start+1) = 2;
    
    m = row_ptr(i+1) + 1;
    color(1,m) = 2;
    color(2,m) = count2;
    color(3,m) = 2;
    color(4,m) = 2;
    
    color(1,m+1) = 2;
    color(2,m+1) = count2;
    color(3,m+1) = 2;
    color(4,m+1) = 1;    
    
    for n = 2:length(fi) 
      % B2 : 2 x (nfe - 2)    
      k = row_start + n;
      color(1,k) = 5; 
      color(2,k) = count2;
      color(3,k) = 1;
      color(4,k) = (n-1);
      
      k = m + n;
      color(1,k) = 5; 
      color(2,k) = count2;
      color(3,k) = 2;
      color(4,k) = (n-1);

      j = fi(n);
      for l = (row_ptr(j) + 1):row_ptr(j+1)
        if col_ind(l) == i         
          % C2 : (nfe - 2) x 2    
          color(1,l) = 8; 
          color(2,l) = count2;
          color(3,l) = n-1;
          color(4,l) = 1;                
          break;
        end
      end
      for l = (row_ptr(j) + 1):row_ptr(j+1)
        if col_ind(l) == (i+1)         
          % C2 : (nfe - 2) x 2    
          color(1,l) = 8; 
          color(2,l) = count2;
          color(3,l) = n-1;
          color(4,l) = 2;                
          break;
        end
      end
    end
    
    idr2(:,count2) = fi(2:end) - nf1;
    idx2(:,:,count2) = nodes2indices(fi(2:end) - nf1, row_ptr2, col_ind2);        
  elseif cnt==3
    % A3 : 3 x 3
    count3 = count3 + 1;
    color(1,row_start) = 3;
    color(2,row_start) = count3;
    color(3,row_start) = 1;
    color(4,row_start) = 1;
    
    color(1,row_start+1) = 3;
    color(2,row_start+1) = count3;
    color(3,row_start+1) = 1;
    color(4,row_start+1) = 2;    
    
    color(1,row_start+2) = 3;
    color(2,row_start+2) = count3;
    color(3,row_start+2) = 1;
    color(4,row_start+2) = 3;    
    
    m = row_ptr(i+1) + 1;
    color(1,m) = 3;
    color(2,m) = count3;
    color(3,m) = 2;
    color(4,m) = 2;
    
    color(1,m+1) = 3;
    color(2,m+1) = count3;
    color(3,m+1) = 2;
    color(4,m+1) = 1;            
    
    color(1,m+2) = 3;
    color(2,m+2) = count3;
    color(3,m+2) = 2;
    color(4,m+2) = 3;            
    
    p = row_ptr(i+2) + 1;
    color(1,p) = 3;
    color(2,p) = count3;
    color(3,p) = 3;
    color(4,p) = 3;
    
    color(1,p+1) = 3;
    color(2,p+1) = count3;
    color(3,p+1) = 3;
    color(4,p+1) = 1;            
    
    color(1,p+2) = 3;
    color(2,p+2) = count3;
    color(3,p+2) = 3;
    color(4,p+2) = 2;        
    
    for n = 3:length(fi) 
      % B3 : 3 x (nfe - 3)    
      k = row_start + n;
      color(1,k) = 6; 
      color(2,k) = count3;
      color(3,k) = 1;
      color(4,k) = (n-2);
      
      k = m + n;
      color(1,k) = 6; 
      color(2,k) = count3;
      color(3,k) = 2;
      color(4,k) = (n-2);
      
      k = p + n;
      color(1,k) = 6; 
      color(2,k) = count3;
      color(3,k) = 3;
      color(4,k) = (n-2);
      
      j = fi(n);
      for l = (row_ptr(j) + 1):row_ptr(j+1)
        if col_ind(l) == i         
          % C3 : (nfe - 3) x 3    
          color(1,l) = 9; 
          color(2,l) = count3;
          color(3,l) = n-2;
          color(4,l) = 1;                
          break;
        end
      end
      for l = (row_ptr(j) + 1):row_ptr(j+1)
        if col_ind(l) == (i+1)         
          % C3 : (nfe - 3) x 3    
          color(1,l) = 9; 
          color(2,l) = count3;
          color(3,l) = n-2;
          color(4,l) = 2;                
          break;
        end
      end
      for l = (row_ptr(j) + 1):row_ptr(j+1)
        if col_ind(l) == (i+2)         
          % C3 : (nfe - 3) x 3    
          color(1,l) = 9; 
          color(2,l) = count3;
          color(3,l) = n-2;
          color(4,l) = 3;                
          break;
        end
      end      
    end
    
    idr3(:,count3) = fi(3:end) - nf1;
    idx3(:,:,count3) = nodes2indices(fi(3:end) - nf1, row_ptr2, col_ind2);        
  end
  
  i = next;
end

idr1 = idr1(:,1:count1);
idr2 = idr2(:,1:count2);
idr3 = idr1(:,1:count3);
idx1 = idx1(:,:,1:count1);
idx2 = idx2(:,:,1:count2);
idx3 = idx1(:,:,1:count3);

for i = (nf1+1):nf % loop over each local face i
  row_start = row_ptr(i) + 1;
  row_end = row_ptr(i+1);      
  color(1,row_start) = 10; % D    
  color(2,row_start) = i;   
  color(3,row_start) = row_ptr2(i-nf1) + 1;   
  color(4,row_start) = row_ptr2(i-nf1) + 1;     
  for k = (row_start+1):row_end
    j = col_ind(k);  % neighboring face j         
    if (j > nf1)  
      color(1,k) = 10; % D    
      color(2,k) = i;         
      for l = (row_ptr2(i-nf1) + 1):row_ptr2(i-nf1+1)
        if col_ind2(l) == j-nf1
          color(3,k) = l;   
          color(4,k) = l;         
        end
      end          
    end    
  end
end

% ind1 = find(f2e(end,:)==0);
% nf1 = length(ind1);
% [count, c] = count_occurrences_sorted(f2e(1,1:nf1));
% if max(count) >=4
%   error("An element on the interface have more than 3 boundary faces.");
% end
% n1 = sum(count==1);
% n2 = sum(count==2);
% n3 = sum(count==3);
% 
% b = zeros(1,nf1);
% for i = 1:nf1
%   if count(i)==1
%     b(i) = i;
%   elseif count(i)==2
%     b(i) = i-n1;
%   elseif count(i)==3
%     b(i) = i-n1-n2;
%   end
% end

% f2e(1,1:nf1)
% count
% b
% c

% nb = row_ptr(end);
% colors = zeros(6,nb);
% for i = 1:nf % loop over each local face i
%   row_start = row_ptr(i) + 1;
%   row_end = row_ptr(i+1);      
%   for k = row_start:row_end
%     j = col_ind(k);  % neighboring face j    
%     colors(5,k) = i;
%     colors(6,k) = j;
%   end
% end
% 
% c1 = zeros(1,nf1);
% c2 = zeros(1,nf1);
% c3 = zeros(1,nf1);
% for i = 1:nf % loop over each local face i
%   row_start = row_ptr(i) + 1;
%   if i <= nf1    
%     colors(1,row_start) = count(i); 
%     % A1 = 1, A2 = 2, A3 = 3
%     if count(i) == 1     % A1 : 1 x 1
%       colors(2,row_start) = i;
%       colors(3,row_start) = i;
%       colors(4,row_start) = i;
%     elseif count(i) == 2 % A2 : 2 x 2
%       colors(2,row_start) = ceil(b(i)/2);
%       colors(3,row_start) = c(i);
%       colors(4,row_start) = c(i);  
%     elseif count(i) == 3 % A3 : 3 x 3
%       colors(2,row_start) = ceil(b(i)/3);
%       colors(3,row_start) = c(i);
%       colors(4,row_start) = c(i);    
%     end
%   else
%     colors(1,row_start) = 10; % D    
%     colors(2,row_start) = i;   
%     colors(3,row_start) = row_ptr2(i-nf1) + 1;   
%     colors(4,row_start) = row_ptr2(i-nf1) + 1;   
%   end
%   b1 = 0;
%   b2 = 0;
%   b3 = 0;
%   row_end = row_ptr(i+1);      
%   for k = (row_start+1):row_end
%     j = col_ind(k);  % neighboring face j         
%     if ( i <= nf1) && (j <= nf1)
%       if colors(1,row_start) == 1 
%         error("something wrong");
%       elseif colors(1,row_start) == 2 % A2 : 2 x 2
%         colors(1,k) = 2;  
%         colors(2,k) = ceil(b(i)/2);
%         colors(3,k) = c(i);
%         colors(4,k) = c(j);          
%       elseif colors(1,row_start) == 3 % A3 : 3 x 3
%         colors(1,k) = 3;    
%         colors(2,k) = ceil(b(i)/3);
%         colors(3,k) = c(i);
%         colors(4,k) = c(j);          
%       end
%     elseif ( i <= nf1) && (j > nf1)
%       % B1 = 4, B2 = 5, B3 = 6
%       if colors(1,row_start) == 1     % B1 : 1 x (nfe-1)
%         b1 = b1 + 1;
%         colors(1,k) = 4;
%         colors(2,k) = i;
%         colors(3,k) = 1;
%         colors(4,k) = b1;
%       elseif colors(1,row_start) == 2 % B2 : 2 x (nfe-2)
%         b2 = b2 + 1;
%         colors(1,k) = 5;  
%         colors(2,k) = ceil(b(i)/2);
%         colors(3,k) = c(i);
%         colors(4,k) = b2;
%       elseif colors(1,row_start) == 3 % B2 : 3 x (nfe-3)
%         b3 = b3 + 1;
%         colors(1,k) = 6;    
%         colors(2,k) = ceil(b(i)/3);
%         colors(3,k) = c(i);
%         colors(4,k) = b3;
%       end            
%     elseif ( i > nf1) && (j <= nf1)  
%       % C1 = 7, C2 = 8, C3 = 9
%       rowj_start = row_ptr(j) + 1;
%       if colors(1,rowj_start) == 1     % C1 : (nfe-1) x 1
%         c1(j) = c1(j) + 1;
%         colors(1,k) = 7;
%         colors(2,k) = j;
%         colors(3,k) = c1(j);
%         colors(4,k) = 1;
%       elseif colors(1,rowj_start) == 2 % C2 : (nfe-2) x 2
%         c2(j) = c2(j) + 1;
%         colors(1,k) = 8;  
%         colors(2,k) = ceil(b(j)/2);
%         colors(3,k) = c2(j);
%         colors(4,k) = c(j);
%       elseif colors(1,rowj_start) == 3 % C3 : (nfe-3) x 3
%         c3(j) = c3(j) + 1;
%         colors(1,k) = 9;    
%         colors(2,k) = ceil(b(j)/3);
%         colors(3,k) = c3(j);
%         colors(4,k) = c(j);
%       end                    
%     elseif ( i > nf1) && (j > nf1)  
%       colors(1,k) = 10; % D    
%       colors(2,k) = i;         
%       for l = (row_ptr2(i-nf1) + 1):row_ptr2(i-nf1+1)
%         if col_ind2(l) == j-nf1
%           colors(3,k) = l;   
%           colors(4,k) = l;         
%         end
%       end      
%     else
%       error("something wrong");
%     end    
%   end
% end

% ind1 = find(f2e(end,:)==0);
% nf1 = length(ind1);
% f2e1 = f2e(:,ind1);
% [row_ptr1, col_ind1] = crs_hdg(f2e1);
% 
% face1 = zeros(nsub,nf1);
% for n = 1:nsub    
%   for i = 1:nf1
%     e = f2e1(1,i);
%     l = f2e1(2,i);
%     face1(n,i) = e2f(l,elem(n,e));
%   end
% end

% ind2 = find(f2e(end,:)>0);
% nf2 = length(ind2);
% f2e2 = f2e(:,ind2);
% [row_ptr2, col_ind2] = crs_hdg(f2e2);

% face2 = zeros(nsub,nf2);
% for n = 1:nsub    
%   for i = 1:nf2
%     e = f2e2(1,i);
%     l = f2e2(2,i);
%     face2(n,i) = e2f(l,elem(n,e));
%   end
% end
% 
end

function [row_ptr, col_ind, f2e, f2f, e2f] = facereorder(elem, t, elemtype, nd)

te = t(:,elem);
[f2e, e2e] = mkf2e(te,elemtype,nd);

nfe = size(e2e,1);
nf = size(f2e,2);

% reorder boundary faces
ind0 = find(f2e(end,:) == 0);
n0 = length(ind0);
f0 = zeros(1,n0);        % count the number of local boundary faces
for i = 1:n0             % loop over each boundary face
  e = f2e(1,ind0(i));    % element associated with this boundary face ind0(i)
  for j = 1:nfe          % loop over local faces of the element e
    if e2e(j,e) == 0     % if there is no neighboring element
      f0(i) = f0(i) + 1; % increment the number of local boundary faces
    end
  end
end
[f0, a] = sort(f0);      % sort the number of local boundary faces
% sort the array a
for j=2:nfe
  k = find(f0==j);       % find faces that have j local boundary faces
  if ~isempty(k)
    % sort elements so that boundary faces on the same element appear consecutively
    [~, b] = sort(f2e(1,ind0(a(k)))); 
    a(k) = a(k(b));
  end
end

ind0 = ind0(a);              % sorted boundary faces
ind1 = setdiff(1:nf, ind0);  % interior faces
f2e = f2e(:,[ind0 ind1]);    % reordered faces

nfsorted = n0;               % count the number of sorted faces
while (nfsorted < nf)        
  e2f = mke2f(f2e);
  [f2f, ~] = mkf2f(f2e, e2f);
  nbf = size(f2f,1);                   % number of neighboring faces  
  unsorted = (nfsorted+1):nf;          % unsorted faces
  nfunsorted = nf - nfsorted;          % number of unsorted faces
  unsortedcount = zeros(1,nfunsorted); % count the number of unsorted neighboring faces  
  for i = 1:nfunsorted                 % loop over each unsorted face
    fi = unsorted(i);                  % unsorted face fi
    for j = 1:nbf                      % loop over each neighboring face of fi
      if f2f(j,fi) > nfsorted          % if the neighboring face is also unsorted
        % increment the number of unsorted neighboring faces
        unsortedcount(i) = unsortedcount(i) + 1;      
      end
    end
  end
  % sort unsortedcount so that unsorted faces with the least unsorted neighboring faces appear first
  [unsortedcount, a] = sort(unsortedcount);    
  % reorder unsorted faces  
  f2e(:,unsorted) = f2e(:,unsorted(a));  
  % update the number of sorted faces
  nfsorted = nfsorted + sum(unsortedcount < nbf); 
end

[row_ptr, col_ind] = crs_hdg(f2e);

% e2f = mke2f(f2e);
% [f2f, f2l] = mkf2f(f2e, e2f);
% 
% row_ptr = zeros(1, nf+1);
% for i = 1:nf
%   n = 1 + sum(f2f(:,i)>0);
%   row_ptr(i+1) = row_ptr(i) + n;
% end
% col_ind = zeros(1, row_ptr(end));
% for i = 1:nf  
%   ind = f2f(:,i)>0;
%   col_ind((row_ptr(i)+1):row_ptr(i+1)) = [i f2f(ind,i)'];
% end

end


% face = zeros(1,nf);
% for i = 1:nf
%   e = f2e(1,i);
%   l = f2e(2,i);
%   face(i) = e2f_global(l,elem(e));
% end

% mesh = mkmesh_square(8,8,1,1,1,1,1,1);
% mesh.xpe=mesh.plocal;
% figure(1); clf; meshplot(mesh,[0 0 1]);
% t = mesh.t;
% elemgroup = 1:16;

