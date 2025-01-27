function [dist,p] = meshdist3(f,dgnodes,perm,ib)

% compute the distance to the wall
perm = perm(:,:,1);
[ir,ic]=find(f==ib);

% compute the nodes on the boundary ib
p = [];
for i = 1:length(ic)
  k = ic(i);
  j = ir(i);  
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
