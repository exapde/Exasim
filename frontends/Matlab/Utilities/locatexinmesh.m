function [elist, xi, y, ind] = locatexinmesh(mesh, master, x, xelem, tol)

porder = mesh.porder;

if size(mesh.t,1) < size(mesh.t,2)    
    mesh.t = mesh.t';
    mesh.p = mesh.p';
end
mesh.nd = master.nd;
mesh.nodetype = master.nodetype;
mesh.elemtype = master.elemtype;

t2e = mkt2e(mesh.t,mesh.elemtype,2);

xdim = ndims(x);
if xdim == 3
  %sz = size(mesh.dgnodes(:,1:mesh.nd,:));  
  sz = size(x);
  x = reshape(permute(x, [1 3 2]), [sz(1)*sz(3) sz(2)]);
end

[elist, xi, y] = findxi(x, xelem, mesh, porder*2, t2e);
ind = findinvalidxi(xi, tol, mesh.elemtype);

if isempty(ind)==0
  for n = 2:6
    [elist2, xi2, y2] = findxi(x(ind,:), elist(ind,:), mesh, porder*(2^n), t2e);
    elist(ind) = elist2;
    xi(ind,:) = xi2;
    y(ind,:) = y2;
    ind = findinvalidxi(xi, tol, mesh.elemtype);
    if isempty(ind)
      break;
    end
  end
end

% xi = permute(reshape(xi, [sz(1) sz(3) sz(2)]),[1 3 2]);
% y = permute(reshape(y, [sz(1) sz(3) sz(2)]), [1 3 2]);
% elist = reshape(elist, [sz(1) sz(3)]);

