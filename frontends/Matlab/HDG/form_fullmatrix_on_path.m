function A = form_fullmatrix_on_path(AE, epath, fpath, fintf, f2e, e2f, elcon)

nep = length(epath);
n1 = size(fintf,1);
n2 = nep*n1 + nep + 1; 

[nfe,ne] = size(e2f);
npf = numel(elcon)/(nfe*ne);
ncu = sqrt(numel(AE)/(ne*npf*npf*nfe*nfe));

elcon = reshape(elcon, [npf nfe ne]);
AE = reshape(AE, [ncu npf nfe ncu npf nfe ne]);

A = zeros(ncu, npf, n2, ncu, npf, n2, 1);

fall = [fintf(:); fpath(1,:)'; fpath(end)]; 
eall = repmat(epath(:)',[n1 1]);
eall = [eall(:); epath(:); epath(end)];
for i = 1:length(fall)
  e = eall(i);
  f = fall(i);
  e1 = f2e(1,f);
  l1 = f2e(2,f);  
  e2 = f2e(3,f);
  l2 = f2e(4,f);            
  j1 = elcon(:,l1,e1) - npf*(f-1);
        
  if e2>0
    j2 = elcon(:,l2,e2) - npf*(f-1);  
    A(:,:,i,:,:,i,1) = AE(:,j1,l1,:,j1,l1,e1) + AE(:,j2,l2,:,j2,l2,e2);          
  else    
    A(:,:,i,:,:,i,1) = AE(:,j1,l1,:,j1,l1,e1);
  end  
  if (i <= nep*n1+1) || (i==n2)    
    if e == e1      
      l12 = l1;      
    elseif e == e2
      l12 = l2;
    end
    j12 = elcon(:,l12,e) - npf*(f-1);       
    for m = 1:nfe
      if m ~= l12
        k = fall == e2f(m,e);
        j3 = elcon(:,m,e) - npf*(e2f(m,e)-1);            
        A(:,:,i,:,:,k,1) = AE(:,j12,l12,:,j3,m,e);
      end
    end    
  else
    for m = 1:nfe
      if m ~= l1
        k = fall == e2f(m,e1);        
        j3 = elcon(:,m,e1) - npf*(e2f(m,e1)-1);    
        A(:,:,i,:,:,k,1) = AE(:,j1,l1,:,j3,m,e1);
      end
      if m ~= l2
        k = fall == e2f(m,e2);
        j3 = elcon(:,m,e2) - npf*(e2f(m,e2)-1);    
        A(:,:,i,:,:,k,1) = AE(:,j2,l2,:,j3,m,e2);
      end
    end    
  end
end

A = reshape(A, ncu*npf*n2, ncu*npf*n2);

% 
% for i = 1:nep  
%   e = epath(i);
%   for j = 1:n1
%     l = lintf(j,i);
%     f = e2f(l,e); 
%     e2 = f2e(3,f);
%     k = j+2*(i-1);
%     if e2>0
%       e1 = f2e(1,f);
%       l1 = f2e(2,f);  
%       l2 = f2e(4,f);      
%       A(:,:,k,:,:,k,1) = AE(:,:,l1,:,:,l1,e1) + AE(:,j2,l2,:,j2,l2,e2);          
%     else
%       A(:,:,k,:,:,k,1) = AE(:,:,l,:,:,l,e);
%     end
%     
%     for q = 1:2
%       m = lpath(q,i);
%       k = nep*n1 ;
%       A(:,:,k,:,:,k,1) = AE(:,:,l,:,:,l,e);
%     end
%   end
% end
% 
% for i = 1:nep  
%   e = epath(i);
%   for j = 1:n1
%     l = lintf(j,i);
%     f = e2f(l,e); 
%     e2 = f2e(3,f);
%     k = j+2*(i-1);
%     if e2>0
%       e1 = f2e(1,f);
%       l1 = f2e(2,f);  
%       l2 = f2e(4,f);      
%       A(:,:,k,:,:,k,1) = AE(:,:,l1,:,:,l1,e1) + AE(:,j2,l2,:,j2,l2,e2);          
%     else
%       A(:,:,k,:,:,k,1) = AE(:,:,l,:,:,l,e);
%     end
%     
%     for q = 1:2
%       m = lpath(q,i);
%       k = nep*n1 ;
%       A(:,:,k,:,:,k,1) = AE(:,:,l,:,:,l,e);
%     end
%   end
% end
%   
