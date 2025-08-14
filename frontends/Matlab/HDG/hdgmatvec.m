function v = hdgmatvec(AE, F, f2e, elcon, opts)

if nargin<5
  opts = 0;
end
  
[npf, nfe, ne] = size(elcon);
ndf = npf*nfe;
nf = size(f2e,2);
n = size(AE,1);
ncu = n/ndf;
ncf = ncu*npf;

Ft = reshape(full(F(:)), [ncu npf*nf]);
FE = reshape(Ft(:,elcon), [ncf*nfe ne]);  
w = zeros(ncf*nfe,ne);
for i = 1:ne  
  w(:,i) = AE(:,:,i)*FE(:,i);  
end
w = reshape(w, [ncu npf nfe ne]);

v = zeros(ncu, npf, nf);
for i = 1:nf % all faces
  e1 = f2e(1,i);
  l1 = f2e(2,i);
  v(:,:,i) = w(:,:,l1,e1);
end

if opts == 0
  for i = 1:nf % interior faces only
    e2 = f2e(3,i);      
    if (e2 > 0)        
      l2 = f2e(4,i);
      ind = elcon(:,l2,e2) - npf*(i-1);    
      v(:,:,i) = v(:,:,i) + w(:,ind,l2,e2);
    end
  end 
end

v = reshape(v, [ncu, npf*nf]);

% for i = 1:ne  
%   for j = 1:nfe
%     vj = w(:,:,j,i);    
%     f = e2f(j, i);    
%     d = elcon(:,j,i) - npf*(f-1);   
%     v(:,:,f) = v(:,:,f) + vj(:,d);
%   end
% end


% AE = permute(reshape(AE,[ncf nfe ncf nfe ne]), [1 3 4 2 5]);
% AE = reshape(AE, [ncf ncf*nfe nfe ne]); 
% 
% Ft = reshape(full(F(:)), [ncu npf*nf]);
% FE = reshape(Ft(:,elcon), [ncf*nfe ne]);  
% F = reshape(FE, [ncf*nfe 1 ne]);
% for j = 1:nfe
%   F(:,j,:) = F(:,1,:);
% end
% 
% w = zeros(ncf,nfe,ne);
% for i = 1:ne  
%   for j = 1:nfe
%     w(:,j,i) = AE(:,:,j,i)*F(:,j,i);
%   end
% end
% w = reshape(w, [ncu npf nfe ne]);

% for i = 1:nf
%   e1 = f2e(1,i);
%   l1 = f2e(2,i);
%   e2 = f2e(3,i);
%     
%   v1 = AE(:,:,l1,e1)*FE(:,e1);  
%   v(:,i) = v(:,i) + v1; 
%   
%   if (e2 > 0)        
%     l2 = f2e(4,i);
%     v2 = reshape(AE(:,:,l2,e2)*FE(:,e2), [ncu npf]);    
%     ind = elcon(:,l2,e2) - npf*(i-1);    
%     v3 = v2(:,ind);
%     v(:,i) = v(:,i) + v3(:);
%   end
% end 
% v = reshape(v, [ncu, npf*nf]);

% F = reshape(FE, [ncf*nfe 1 ne]);
% for j = 1:nfe
%   F(:,j,:) = F(:,1,:);
% end
% w = zeros(ncf,nfe,ne);
% for i = 1:ne  
%   for j = 1:nfe
%     w(:,j,i) = AE(:,:,j,i)*F(:,j,i);
%   end
% end
% w = reshape(w, [ncu npf nfe ne]);
% v = zeros(ncu, npf, nf);
% for i = 1:nf
%   e1 = f2e(1,i);
%   l1 = f2e(2,i);
%   v(:,:,i) = v(:,:,i) + w(:,:,l1,e1);
%   e2 = f2e(3,i);      
%   if (e2 > 0)        
%     l2 = f2e(4,i);
%     ind = elcon(:,l2,e2) - npf*(i-1);    
%     v(:,:,i) = v(:,:,i) + w(:,ind,l2,e2);
%   end
% end 
% for i = 1:ne  
%   for j = 1:nfe
%     vj = w(:,:,j,i);    
%     f = e2f(j, i);    
%     d = elcon(:,j,i) - npf*(f-1);   
%     v(:,:,f) = v(:,:,f) + vj(:,d);
%   end
% end


% for i = 1:ne  
%   for j = 1:nfe
%     vj = AE(:,:,j,i)*FE(:,i);    
%     f = e2f(j, i);    
%     d = elcon(:,j,i) - npf*(f-1);   
%     v(:,f) = v(:,f) + vj(d);
%   end
% end
% v = reshape(v, [ncu, npf*nf]);



% Ft = reshape(full(F(:)), [ncu npf*nf]);
% FE = reshape(Ft(:,elcon), [ncf nfe ne]);  
% AE  = permute(reshape(AE,[ncf nfe ncf nfe ne]), [1 3 2 4 5]);
% v = zeros(ncf, nf);
% 
% for i = 1:nf
%   e1 = f2e(1,i);
%   l1 = f2e(2,i);
%   e2 = f2e(3,i);
%   
%   v1 = AE(:,:,l1,l1,e1)*FE(:,l1,e1);
%   for j = 1:nfe
%     if j ~= l1
%       v1 = v1 + AE(:,:,l1,j,e1)*FE(:,j,e1);
%     end
%   end    
%   v(:,i) = v(:,i) + v1; 
%   
%   if (e2 > 0)    
%     l2 = f2e(4,i);
%     v2 = AE(:,:,l2,l2,e2)*FE(:,l2,e2);
%     for j = 1:nfe
%       if j ~= l2
%         v2 = v2 + AE(:,:,l2,j,e2)*FE(:,j,e2);
%       end
%     end
%     ind = elcon(:,l2,e2) - npf*(i-1);
%     v(:,i) = v(:,i) + v2(ind);
%   end
% end
%  
% v = reshape(v, [ncu, npf*nf]);


% for i = 1:ne
%   fe = FE(:,:,i);
%   ae = AE(:,:,:,:,i);    
%   for j = 1:nfe
%     vj = zeros(ncf,1);
%     for k = 1:nfe      
%       vj = vj + ae(:,:,j,k)*fe(:,k);
%     end
%     f = e2f(j, i);
%     d = elcon(:,j,i) - npf*(f-1);
%     v(:,f) = v(:,f) + vj(d);
%   end
% end

