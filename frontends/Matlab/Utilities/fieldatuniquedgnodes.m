function UDG1 = fieldatuniquedgnodes(mesh, master, UDG, dgnodes1)

[npe, ~, ne]  = size(dgnodes1);
nc = size(UDG, 2);
nd = size(dgnodes1,2);

if nd == 1
  x = dgnodes1(:,1,:); 
  p = [x(:)];    
elseif nd == 2
  x = dgnodes1(:,1,:);
  y = dgnodes1(:,2,:);  
  p = [x(:) y(:)];    
elseif nd == 3
  x = dgnodes1(:,1,:);
  y = dgnodes1(:,2,:);
  z = dgnodes1(:,2,:);
  p = [x(:) y(:) z(:)];    
end

[pu, ~, ip] = unique(p, 'rows', 'stable');

[elist, xi] = locatexinmesh(mesh, master, pu, [], 1e-4);
udgu = evalfield(mesh, UDG, elist, xi);

UDG1 = udgu(ip,:);
UDG1 = permute(reshape(UDG1, [npe ne nc]),[1 3 2]); 

