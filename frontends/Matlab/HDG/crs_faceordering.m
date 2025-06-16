function [row_ptr, col_ind, face, f2eelem] = crs_faceordering(elem, f2e)

ne1 = max(f2e(1,:));
ne2 = max(f2e(3,:));
ne = max(ne1, ne2); 
nfe1 = max(f2e(2,:));
nfe2 = max(f2e(4,:));
nfe = max(nfe1, nfe2); 

e2f = mke2f(f2e, nfe, ne);
nfe = size(e2f,1);
[nb, ne] = size(elem);

[felem, nfelem] = mkfelem(e2f, elem(1,:), nfe, ne);

f2eelem = mkf2eelem(f2e, felem, elem(1,:), nfelem, ne);

f2eelem = facereorder(f2eelem);

[row_ptr, col_ind] = crs_array(f2eelem, nfe, ne, nfelem);

face = mkface(e2f, elem, f2eelem, nb, nfelem);

% e2felem = e2f(:,elem(1,:));
% [nfe, ne] = size(e2felem);
% felem = zeros(1,nfe*ne);
% n = 0;
% for i = 1:ne
%   for j = 1:nfe
%     f = e2felem(j,i);
%     found = 0;
%     for k = 1:n
%       if felem(k) == f
%         found = 1;
%         break;
%       end
%     end
%     if found == 0
%       n = n + 1;
%       felem(n) = f;
%     end
%   end  
% end

% f2eelem = f2e(:,felem(1:n));
% for i = 1:n
%   e1 = f2eelem(1,i);
%   e2 = f2eelem(3,i);
%   
%   found1 = 0;
%   for j = 1:ne
%     if elem(1,j) == e1    
%       found1 = 1;
%       f2eelem(1,i) = j;
%       break;
%     end
%   end
%   
%   found2 = 0;
%   for j = 1:ne
%     if elem(1,j) == e2    
%       found2 = 1;
%       f2eelem(3,i) = j;
%       break;
%     end
%   end
%   
%   if (found1==0) && (found2==1)
%     f2eelem(1,i) = f2eelem(3,i);
%     f2eelem(2,i) = f2eelem(4,i);
%     f2eelem(3,i) = 0;
%     f2eelem(4,i) = 0;
%   end  
%   if found2==0
%     f2eelem(3,i) = 0;
%     f2eelem(4,i) = 0;
%   end  
% end

% nb = size(elem,1);
% nf = size(f2eelem,2);
% face = zeros(nb,nf);
% for n = 1:nb    
%   for i = 1:nf
%     e = f2eelem(1,i);
%     l = f2eelem(2,i);
%     face(n,i) = e2f(l,elem(n,e));
%   end
% end

end

function [felem,n] = mkfelem(e2f, elem, nfe, ne)

e2felem = e2f(:,elem);
felem = zeros(1,nfe*ne);
n = 0;
for i = 1:ne
  for j = 1:nfe
    f = e2felem(j,i);
    found = 0;
    for k = 1:n
      if felem(k) == f
        found = 1;
        break;
      end
    end
    if found == 0
      n = n + 1;
      felem(n) = f;
    end
  end  
end

end

function f2eelem = mkf2eelem(f2e, felem, elem, n, ne)

f2eelem = f2e(:,felem(1:n));
for i = 1:n
  e1 = f2eelem(1,i);
  e2 = f2eelem(3,i);
  
  found1 = 0;
  for j = 1:ne
    if elem(j) == e1    
      found1 = 1;
      f2eelem(1,i) = j;
      break;
    end
  end
  
  found2 = 0;
  for j = 1:ne
    if elem(j) == e2    
      found2 = 1;
      f2eelem(3,i) = j;
      break;
    end
  end
  
  if (found1==0) && (found2==1)
    f2eelem(1,i) = f2eelem(3,i);
    f2eelem(2,i) = f2eelem(4,i);
    f2eelem(3,i) = 0;
    f2eelem(4,i) = 0;
  end  
  if found2==0
    f2eelem(3,i) = 0;
    f2eelem(4,i) = 0;
  end  
end

end


function f2e = facereorder(f2e)

ind0 = find(f2e(end,:) == 0);
[~, ii] = sortrecurrence(f2e(1,ind0));
ind0 = ind0(ii);

nf = size(f2e,2);
ind1 = setdiff(1:nf, ind0);  % interior faces

f2e = f2e(:,[ind0 ind1]);    % reordered faces

ne1 = max(f2e(1,:));
ne2 = max(f2e(3,:));
ne = max(ne1, ne2); 
nfe1 = max(f2e(2,:));
nfe2 = max(f2e(4,:));
nfe = max(nfe1, nfe2); 

f2e = sortinteriorfaces(f2e, nfe, ne, nf, length(ind0));

end

% nf = size(f2e,2);
% e2f = mke2f(f2e);
% e2e = mke2e(f2e, e2f);
% nfe = size(e2f,1);
% % reorder boundary faces
% ind0 = find(f2e(end,:) == 0);
% n0 = length(ind0);
% f0 = zeros(1,n0);        % count the number of local boundary faces
% for i = 1:n0             % loop over each boundary face
%   e = f2e(1,ind0(i));    % element associated with this boundary face ind0(i)
%   for j = 1:nfe          % loop over local faces of the element e
%     if e2e(j,e) == 0     % if there is no neighboring element
%       f0(i) = f0(i) + 1; % increment the number of local boundary faces
%     end
%   end
% end
% 
% [f0, a] = sort(f0);      % sort the number of local boundary faces
% % sort the array a
% for j=2:nfe
%   k = find(f0==j);       % find faces that have j local boundary faces
%   if ~isempty(k)
%     % sort elements so that boundary faces on the same element appear consecutively
%     [~, b] = sort(f2e(1,ind0(a(k)))); 
%     a(k) = a(k(b));
%   end
% end
% 
% ind0 = ind0(a);              % sorted boundary faces

% ind0 = find(f2e(end,:) == 0);
% [~, ii] = sortrecurrence(f2e(1,ind0));
% ind0 = ind0(ii);
% 
% nf = size(f2e,2);
% ind1 = setdiff(1:nf, ind0);  % interior faces
% f2e = f2e(:,[ind0 ind1]);    % reordered faces
% 
% ne1 = max(f2e(1,:));
% ne2 = max(f2e(3,:));
% ne = max(ne1, ne2); 
% nfe1 = max(f2e(2,:));
% nfe2 = max(f2e(4,:));
% nfe = max(nfe1, nfe2); 
% 
% f2e = sortinteriorfaces(f2e, nfe, ne, nf, length(ind0));

% nfsorted = length(ind0);               % count the number of sorted faces
% while (nfsorted < nf)        
%   e2f = mke2f(f2e, nfe, ne);
%   f2f = mkf2f(f2e, e2f, nfe, nf);
%   nbf = 2*(nfe-1);                     % maximum number of neighboring faces  
%   unsorted = (nfsorted+1):nf;          % unsorted faces
%   nfunsorted = nf - nfsorted;          % number of unsorted faces
%   unsortedcount = zeros(1,nfunsorted); % count the number of unsorted neighboring faces  
%   for i = 1:nfunsorted                 % loop over each unsorted face
%     fi = unsorted(i);                  % unsorted face fi
%     for j = 1:nbf                      % loop over each neighboring face of fi
%       if f2f(j,fi) > nfsorted          % if the neighboring face is also unsorted
%         % increment the number of unsorted neighboring faces
%         unsortedcount(i) = unsortedcount(i) + 1;      
%       end
%     end
%   end
%   
%   % sort unsortedcount so that unsorted faces with the least unsorted neighboring faces appear first
%   [unsortedcount, a] = simple_bubble_sort(unsortedcount);    
% 
%   % reorder unsorted faces  
%   f2e(:,unsorted) = f2e(:,unsorted(a));  
%   % update the number of sorted faces
%   nfsorted = nfsorted + sum(unsortedcount < nbf); 
% end

% end

function f2e = sortinteriorfaces(f2e, nfe, ne, nf, nf0)

nfsorted = nf0;                        % count the number of sorted faces
while (nfsorted < nf)        
  e2f = mke2f(f2e, nfe, ne);
  f2f = mkf2f(f2e, e2f, nfe, nf);
  nbf = 2*(nfe-1);                     % maximum number of neighboring faces  
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
  [unsortedcount, a] = simple_bubble_sort(unsortedcount);    

  % reorder unsorted faces  
  f2e(:,unsorted) = f2e(:,unsorted(a));  
  % update the number of sorted faces
  nfsorted = nfsorted + sum(unsortedcount < nbf);     
end

end

function face = mkface(e2f, elem, f2eelem, nb, nf)

face = zeros(nb,nf);
for n = 1:nb    
  for i = 1:nf
    e = f2eelem(1,i);
    l = f2eelem(2,i);
    face(n,i) = e2f(l,elem(n,e));
  end
end

end

function e2f = mke2f(f2e, nfe, ne)

e2f = zeros(nfe, ne);
nf = size(f2e, 2);
for i = 1:nf
  e1 = f2e(1,i);
  l1 = f2e(2,i);
  e2 = f2e(3,i);
  e2f(l1, e1) = i;  
  if e2>0
    l2 = f2e(4,i);
    e2f(l2, e2) = i;
  end
end
end

function [f2f, f2l] = mkf2f(f2e, e2f, nfe, nf)
%  
%   f2e   :  Face to element connectivity
%   e2f   :  Element to face connectivity
%
%   f2f   :  Face to face connectivity

nbf = 2*(nfe-2);    % number of neighboring faces
f2f = zeros(nbf,nf);
f2l = zeros(nbf,nf);
for i = 1:nf  % loop over each face
    e1 = f2e(1,i);
    l1 = f2e(2,i);
    e2 = f2e(3,i);
    l2 = f2e(4,i);    
    
    k = 0;       
    for l=1:nfe % loop over each faces of the 1st element
        if l ~= l1  
            k = k + 1;                 
            j = e2f(l,e1);
            f2f(k,i) = j;
            f2l(k,i) = l;
        end
    end
    
    if e2>0           % face i is an interior face                
        for l=1:nfe % loop over faces of the 2nd element
            if l ~= l2                                                
                k = k + 1;                 
                j = e2f(l,e2);
                f2f(k,i) = j;
                f2l(k,i) = l;
            end
        end        
    end    
end

end

function [row_ptr, col_ind] = crs_array(f2e, nfe, ne, nf)  

  e2f = mke2f(f2e, nfe, ne);
  f2f = mkf2f(f2e, e2f, nfe, nf);

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

end

function [b, e] = sortrecurrence(a)

n = length(a);
recurrence = zeros(1, n);

% Step 1: Count occurrences for each entry
for i = 1:n
    count = 0;
    for j = 1:n
        if a(i) == a(j)
            count = count + 1;
        end
    end
    recurrence(i) = count;
end

maxrec = 0;
for i = 1:n
    if recurrence(i) > maxrec
        maxrec = recurrence(i);
    end
end

b = zeros(1, n);
e = zeros(1, n);
c = zeros(1, n);
d = zeros(1, n);
m = 0;
for k = 1:maxrec
    j = 0;    
    for i = 1:n
        if recurrence(i) == k
            j = j + 1;
            c(j) = a(i);        
            d(j) = i;
        end
    end
    [c, ind] = simple_bubble_sort(c(1:j));
    d = d(ind);
    b(m+1:m+j) = c;
    e(m+1:m+j) = d;
    m = m + j;
    if m>=n
      break;
    end
end

if max(abs(b-a(e))) > 0
  error("something wrong");
end

end

function [b, ind] = simple_bubble_sort(a)
n = length(a);
b = a;
ind = 1:n;

for i = 1:n-1
    for j = 1:n-i
        if b(j) > b(j+1)
            % Swap values
            tmp = b(j);
            b(j) = b(j+1);
            b(j+1) = tmp;
            % Swap indices
            tmp_idx = ind(j);
            ind(j) = ind(j+1);
            ind(j+1) = tmp_idx;
        end
    end
end
end

