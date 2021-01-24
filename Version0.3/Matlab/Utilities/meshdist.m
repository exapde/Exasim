function dist = meshdist(f,dgnodes,perm,ib)
% compute the distance to the wall
perm = perm(:,:,1);
ind = [];
for i = 1:length(ib)
    ind = [ind; find(f(end,:)==-ib(i))];
end
% compute the nodes on the boundary ib
p = [];
for ii=ind
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

% function dist = meshdist(mesh,ib)
% % compute the distance to the wall
% perm = mesh.perm(:,:,1);
% ind = [];
% for i = 1:length(ib)
%     ind = [ind; find(mesh.f(:,end)==-ib(i))];
% end
% % compute the nodes on the boundary ib
% nfe = size(perm,2);
% p = [];
% for ii=ind'
%   k = mesh.f(ii,end-1);
%   %fc = mesh.bf(:,k);
%   %bf = (fc<0);
%   for j=1:nfe
%       %if bf(j)==1
%       if mesh.t2f(k,j)==ii
%         p = [p; mesh.dgnodes(perm(:,j),:,k)];
%       end
%   end
% end
% nn = size(mesh.dgnodes,1);
% ne = size(mesh.dgnodes,3);
% dist = zeros(nn,1,ne);
% for i = 1:ne
%     for j=1:nn
%         x = mesh.dgnodes(j,:,i);
%         s = sqrt((p(:,1)-x(1)).^2 + (p(:,2)-x(2)).^2);
%         dist(j,1,i) = min(s);
%     end
% end
