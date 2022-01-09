function [p,t] = removeelemementinpolygon(p, t, xv, yv)

ne = size(t,1); % # elements
nv = size(t,2); % # vertices
nd = size(p,2); % # dimensions

p1 = reshape(p(t',:),[nv ne nd]);
p1 = reshape(mean(p1,1),[ne nd]);

ind = inpolygon(p1(:,1),p1(:,2),xv,yv);
t(ind,:) = [];

[p,t]=fixmesh(p,t);




