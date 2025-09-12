function [xi, y] = newtonx(dgnodes, x, elist, xilist, plocal, elemtype, porder)

npe = size(dgnodes, 1);
nd = size(dgnodes, 2);
ne = size(xilist,1);

y = zeros(ne, nd);
dydxi = zeros(ne, nd, nd);
A = zeros(ne, nd, nd);

dgnodes = permute(dgnodes(:,:,elist),[3 1 2]); % ne x npe x nd
xi = 1.0*xilist;

for iter = 1:5
  nfs = mkshape(porder, plocal, xi, elemtype); % npe x ne x (nd+1)
  nfs = permute(nfs, [2 1 3]);        % ne x npe x (nd+1)
  y(:,1) = dgnodes(:,1,1).*nfs(:,1,1);  % ne x nd
  y(:,2) = dgnodes(:,1,2).*nfs(:,1,1);  % ne x nd
  dydxi(:,1,1) = dgnodes(:,1,1).*nfs(:,1,2);
  dydxi(:,1,2) = dgnodes(:,1,1).*nfs(:,1,3);
  dydxi(:,2,1) = dgnodes(:,1,2).*nfs(:,1,2);
  dydxi(:,2,2) = dgnodes(:,1,2).*nfs(:,1,3);
  for n = 2:npe    
    y(:,1) = y(:,1) + dgnodes(:,n,1).*nfs(:,n,1); 
    y(:,2) = y(:,2) + dgnodes(:,n,2).*nfs(:,n,1); 
    dydxi(:,1,1) = dydxi(:,1,1) + dgnodes(:,n,1).*nfs(:,n,2);
    dydxi(:,1,2) = dydxi(:,1,2) + dgnodes(:,n,1).*nfs(:,n,3);
    dydxi(:,2,1) = dydxi(:,2,1) + dgnodes(:,n,2).*nfs(:,n,2);
    dydxi(:,2,2) = dydxi(:,2,2) + dgnodes(:,n,2).*nfs(:,n,3);
  end

%   if iter == 1
%     e = y - xlist;
%     if max(e(:,1).^2 + e(:,2).^2) > 1e-12
%       error("something wrong");
%     end
%   end
  
  r = x - y;  
  detJ = dydxi(:,1,1).*dydxi(:,2,2) - dydxi(:,1,2).*dydxi(:,2,1);  
  A(:,1,1) = dydxi(:,2,2)./detJ;
  A(:,1,2) = -dydxi(:,1,2)./detJ;
  A(:,2,1) = -dydxi(:,2,1)./detJ;
  A(:,2,2) = dydxi(:,1,1)./detJ; 
  dxi(:,1) = A(:,1,1).*r(:,1) + A(:,1,2).*r(:,2);
  dxi(:,2) = A(:,2,1).*r(:,1) + A(:,2,2).*r(:,2);
  xi = xi + dxi;    
end


