function [elist, xi, y] = findxi(x, xelem, mesh, pgauss, t2e)
% x : input points 
% xlem : element index for each point in x 
% mesh : the high-order mesh
% pgauss : the degree for Gauss points
% t2e : element-to-element connectivities 

% elist : elist(n) is the element containing the point x(n, :) 
% xi : xi(n,:) is the master point corresponding to the point x(n, :) 
% y : y(n,:) is the physical point corresponding to the point x(n, :) 
% y must be equal to x

master = mkmaster(mesh, pgauss);
 
% reshape x if necessary
xdim = ndims(x);
if xdim == 3
  szx = size(x);
  x = reshape(permute(x,[1 3 2]), [szx(1)*szx(3) szx(2)]);
  if isempty(xelem)==1
    xelem = repmat((1:szx(3)),[szx(1) 1]);
  end  
end

if isempty(xelem)==1
  % element centers
  nd = mesh.nd;  
  [ne,nv] = size(mesh.t);
  p = reshape(mesh.p(mesh.t',:),[nv ne nd]);
  xm = reshape(mean(p,1),[ne nd]);
  nx = size(x,1);
  xelem = zeros(nx,1);  
  for i = 1:nx
    d = (x(i,1)-xm(:,1)).^2 + (x(i,2)-xm(:,2)).^2;        
    [~,xelem(i)] = min(d);
  end
end

% For each point in x, find its nearest gauss point and the element
% containing the gauss points
[elist, ~, xilist] = findex(mesh.dgnodes(:,1:mesh.nd,:), t2e, x, xelem, master.gpvl, master.shapmv);

% use Newton method to find xi for each point in x
[xi, y] = newtonx(mesh.dgnodes(:,1:mesh.nd,:), x, elist, xilist, master.plocvl, mesh.elemtype, master.porder);

if xdim == 3
  xi = permute(reshape(xi, [szx(1) szx(3) szx(2)]),[1 3 2]);
  y = permute(reshape(y, [szx(1) szx(3) szx(2)]), [1 3 2]);
  elist = reshape(elist, [szx(1) szx(3)]);
end


% % find all invalid points
% tol = 1e-4;
% ind = findind(xi, tol, mesh.elemtype);
% 

% xg = master.shapmv(:,:,1)*reshape(mesh.dgnodes, master.npv,[]);
% xg = reshape(xg, [master.ngv master.nd mesh.ne]);
% ex = elist(ind);
% e = t2e(ex,:);
% k = min(e(:,1))-1;
% for j = 1:k  
%   ex = elist(ind);
%   e = t2e(ex,:);  
%   ej = e(:,j+2);
%   n = length(ind);
%   for i = 1:n
%     s = (x(ind(i),1) - xg(:,1,ej(i))).^2 + (x(ind(i),2) - xg(:,2,ej(i))).^2;
%     [~,ii] = min(s);
%     xlist(i,:) = xg(ii,:,ej(i));
%     xilist(i,:) = master.gpvl(ii,:);       
%   end
%   [a, b] = newtonx(mesh.dgnodes, x(ind,:), ej, xlist(1:n,:), xilist(1:n,:), master.plocvl, mesh.elemtype, master.porder);
%   
%   ind2 = findind(a, tol, mesh.elemtype);      
%   ind1 = setdiff(1:size(a,1), ind2);     
%   
%   xi(ind(ind1),:) = a(ind1,:);
%   y(ind(ind1),:) = b(ind1,:);
%   elist(ind(ind1)) = ej(ind1);              
%   ind = ind(ind2);    
% end
% 
% ind = findind(xi, tol, mesh.elemtype);
% for i = 1:length(ind)
%   z = x(ind(i),:);
%   e = elist(ind(i));
%   e = t2e(e,2:end);
%   e = e(e>0);
%   e = t2e(e,2:end);
%   e = unique(e(:));
%   e = e(e>0);  
%   
%   n = length(e);  
%   s = (z(1) - xg(:,1,e)).^2 + (z(2) - xg(:,2,e)).^2;
%   smin = min(s(:));
%   s = s - smin;
%   [ii,jj] = find(s==0);
%   ei = e(jj);  
%   zi = repmat(xg(ii,:,ei),[n 1]);
%   vi = repmat(master.gpvl(ii,:), [n 1]);
%     
%   [ai, bi] = newtonx(mesh.dgnodes, repmat(z, [n 1]), e(:), zi, vi, master.plocvl, mesh.elemtype, master.porder);
%   ind2 = findind(ai, tol, mesh.elemtype);  
%   ind1 = setdiff(1:size(ai,1), ind2);     
%   ind1 = ind1(1);
%   xi(ind(i),:) = ai(ind1,:);
%   y(ind(i),:) = bi(ind1,:);
%   elist(ind(i)) = e(ind1);            
% end

% ind = findind(xi, tol, mesh.elemtype);
% length(ind)


