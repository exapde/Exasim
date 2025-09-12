function [A, B, C, D] = pathcompute(AE, epath, fpath, lpath, fintf, lintf, nelems, f2e, e2f, elcon)

[A, B, C, D] = pathsystem(AE, epath, fpath, lpath, fintf, lintf, f2e, e2f, elcon);

npaths = length(nelems)-1;
for n = 1:npaths
  ind = (nelems(n)+1):nelems(n+1);
  [A(:,:,ind), C(:,:,ind), D(:,:,ind)] = pathlu(A(:,:,ind), B(:,:,ind), C(:,:,ind), D(:,:,ind));    
end

% ne = size(A,3);
% for i = 1:ne
%     A(:,:,i) = inv(A(:,:,i));    
%     C(:,:,i) = C(:,:,i)*A(:,:,i);    
% end
% 
% m = size(D, 1)/2;
% 
% for n = 1:npaths
%   ind = (nelems(n)+1):nelems(n+1);
%   nep = length(ind);
%   for k = 1:nep
%     i = ind(k);
%     CB = C(:,:,i)*B(:,:,i);
%     D(:,:,i) = D(:,:,i) - CB;
%     if k < nep    
%       D(1:m,1:m,i+1) = D(1:m,1:m,i+1) - CB(m+1:2*m,m+1:2*m); 
%     end
%   end
% end
% 
% nep = nelems(2)-nelems(1);
% % LU decomposition of tridagonal system
% for k = 1:nep
%   for n = 1:npaths
%     ind = (nelems(n)+1):nelems(n+1);
%     i = ind(k);
%     D(1:m,1:m,i) = inv(D(1:m,1:m,i));
%     D(m+1:2*m,1:m,i) = D(m+1:2*m,1:m,i)*D(1:m,1:m,i);   
%     if k < nep
%       D(1:m,1:m,i+1) = D(1:m,1:m,i+1) - D(m+1:2*m,1:m,i)*D(1:m,m+1:2*m,i);
%     else
%       D(m+1:2*m,m+1:2*m,i) = D(m+1:2*m,m+1:2*m,i) - D(m+1:2*m,1:m,i)*D(1:m,m+1:2*m,i);
%       D(m+1:2*m,m+1:2*m,i) = inv(D(m+1:2*m,m+1:2*m,i));
%     end
%   end  
% end
% 
% % for n = 1:npaths
% %   ind = (nelems(n)+1):nelems(n+1);
% %   nep = length(ind);
% %   
% %   % LU decomposition of tridagonal system
% %   for k = 1:nep
% %     i = ind(k);
% %     D(1:m,1:m,i) = inv(D(1:m,1:m,i));
% %     D(m+1:2*m,1:m,i) = D(m+1:2*m,1:m,i)*D(1:m,1:m,i);   
% %     if k < nep
% %       D(1:m,1:m,i+1) = D(1:m,1:m,i+1) - D(m+1:2*m,1:m,i)*D(1:m,m+1:2*m,i);
% %     else
% %       D(m+1:2*m,m+1:2*m,i) = D(m+1:2*m,m+1:2*m,i) - D(m+1:2*m,1:m,i)*D(1:m,m+1:2*m,i);
% %       D(m+1:2*m,m+1:2*m,i) = inv(D(m+1:2*m,m+1:2*m,i));
% %     end
% %   end  
% % end
% 
% 
