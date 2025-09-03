function uxi = evalfield(master, u, elem, xi)

%u = permute(u,[3 1 2]); % ne x npe x m

nfs = mkshape(master.porder, master.xpe, xi, master.elemtype); % npe x n x (nd+1)
nfs = permute(nfs(:,:,1), [2 1]); % n x npe 

n = size(xi,1);
m = size(u,2);

uxi = zeros(n,m);
for i = 1:n
  e = elem(i);  
  uxi(i,:) = nfs(i,:)*u(:,:,e); 
end


% M = size(u,3);
% N = size(xi,1);
% y = zeros(N,M);
% 
% for m = 1:M
%   y(:,m) = u(:,1,m).*nfs(:,1,1);  % ne x nd
% end
% 
% y(:,2) = dgnodes(:,1,2).*nfs(:,1,1);  % ne x nd  
% for n = 2:npe    
%   y(:,1) = y(:,1) + dgnodes(:,n,1).*nfs(:,n,1); 
%   y(:,2) = y(:,2) + dgnodes(:,n,2).*nfs(:,n,1); 
% end
  


