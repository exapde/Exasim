function [p,t] = rmelems(p, t, pc, rc)

ne = size(t,1); % # elements
nv = size(t,2); % # vertices
nd = size(p,2); % # dimensions

p1 = reshape(p(t',:),[nv ne nd]);
p1 = reshape(mean(p1,1),[ne nd]);


d = (p1(:,1)-pc(1)).^2;
for j = 2:nd
    d = d + (p1(:,j)-pc(j)).^2;
end
d = sqrt(d);

ind = d<=rc;
t(ind,:) = [];

[p,t]=fixmesh(p,t);




