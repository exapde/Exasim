function [p,t] = removeelemement(p, t, expr)

ne = size(t,1); % # elements
nv = size(t,2); % # vertices
nd = size(p,2); % # dimensions

p1 = reshape(p(t',:),[nv ne nd]);
p1 = reshape(mean(p1,1),[ne nd]);

if nd==1
    x = p1(:,1);    
elseif nd==2
    x = p1(:,1);
    y = p1(:,2);
else
    x = p1(:,1);
    y = p1(:,2);
    z = p1(:,3);
end
   
ind = eval(expr);
t(ind,:) = [];

[p,t]=fixmesh(p,t);




