function [elist, xlist, xilist] = findex(dgnodes, t2e, x, elist, xig, shapmv)
% For each point in x, find its nearest gauss point and the element
% containing the gauss points

nx = size(x,1);
[npe, nd, ne] = size(dgnodes);
nge = size(shapmv, 1);

% Gauss points on the entire mesh
xg = shapmv(:,:,1)*reshape(dgnodes, npe,[]);
xg = reshape(xg, [nge nd ne]);

xlist = zeros(nx,nd);
xilist = zeros(nx,nd);
for i = 1:nx
  smin = 0; sm = 1;
  iter = 0;
  while (smin<sm) 
    iter = iter + 1;
    if iter>1
      sm = smin;
    end
    e = elist(i);         % the element associated with x(i,:)
    k = t2e(e, 1);
    el = t2e(e, 2:(k+1)); % list of elements around e
    p = xg(:,:,el);       % Gauss points on those elements 
    s = (x(i,1) - p(:,1,:)).^2 + (x(i,2) - p(:,2,:)).^2;
    if nd == 3
      s = s + (x(i,3) - p(:,3,:)).^2;
    end
    s = reshape(s,[nge k]);
    smin = min(s(:));
    s = s - smin;    
    [a,b] = find(s==0, 1);      
    elist(i) = el(b);
    xlist(i,:) = p(a,:,b);
    xilist(i,:) = xig(a,:);
  end  
end

% x = rand(5,2)-0.5
% elist = zeros(size(x,1),1);
% xg = master.shapmv(:,:,1)*reshape(mesh.dgnodes,master.npv,[]);
% xg = reshape(xg, [master.ngv master.nd mesh.ne]);
% for i = 1:size(x,1)
%   s = (x(i,1) - xg(:,1,:)).^2 + (x(i,2) - xg(:,2,:)).^2;
%   s = squeeze(s);
%   s = s - min(s(:));
%   [a,b] = find(s==0);  
%   elist(i) = b;
% end

