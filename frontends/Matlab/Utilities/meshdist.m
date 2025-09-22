function dist = meshdist(f,dgnodes,perm,ib)

if nargin == 2
  dist = meshdist2(f,dgnodes);
  return;
end

% compute the distance to the wall
perm = perm(:,:,1);
ind = [];
for i = 1:length(ib)
    tm = find(f(end,:)==-ib(i));
    ind = [ind; tm(:)];
end

% compute the nodes on the boundary ib
p = [];
for i = 1:length(ind)
  ii = ind(i);  
  k = f(end-3,ii);
  j = f(end-1,ii);  
  p = [p; dgnodes(perm(:,j),:,k)];  
end
nn = size(dgnodes,1);
ne = size(dgnodes,3);
dist = zeros(nn,1,ne);
for i = 1:ne
    for j=1:nn
        x = dgnodes(j,:,i);
        s = sqrt((p(:,1)-x(1)).^2 + (p(:,2)-x(2)).^2);
        dist(j,1,i) = min(s);
    end
end

end

function dist = meshdist2(mesh,ib)
% compute the distance to the wall

nd = mesh.nd;
perm = mesh.perm(:,:,1);
ind  = find(mesh.f(:,end)==-ib(1));
for l = 2:length(ib)
  tm  = find(mesh.f(:,end)==-ib(l));
  ind = [ind; tm];
end

% compute the nodes on the boundary ib
nfe = size(perm,2);
p = [];
for ii=ind'
  k = mesh.f(ii,end-1);  
  fc = mesh.bf(:,k);
  %bf = (fc<0);  
  for j=1:nfe
      if ismember(fc(j), -ib) %fc(j)==-ib(1)
        p = [p; mesh.dgnodes(perm(:,j),1:nd,k)];        
      end
  end
end

nn = size(mesh.dgnodes,1);
ne = size(mesh.dgnodes,3);
dist = zeros(nn,1,ne);
for i = 1:ne
    for j=1:nn
        x = mesh.dgnodes(j,1:nd,i);
        s = sqrt((p(:,1)-x(1)).^2 + (p(:,2)-x(2)).^2);
        dist(j,1,i) = min(s);  
    end
end

end
