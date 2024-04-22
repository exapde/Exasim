function [p,t] = mergeelemementatnode(p, t, pv)

ne = size(t,1); % # elements
nv = size(t,2); % # vertices
nd = size(p,2); % # dimensions

ind2 = zeros(size(pv,1),1);
for i = 1:size(pv,1)
    d = (p(:,1) - pv(i,1)).^2;
    for j = 2:nd
        d = d + (p(:,j) - pv(i,j)).^2;
    end
    [~,idx] = min(d);
    %[p(idx,:) pv(i,:)]
    [ind,~] = find(t==idx);        
    if length(ind)==2
        ii = find(t(ind(1),:)==idx);
        for jj = 1:nv
            if min(abs(t(ind(2),jj)-t(ind(1),:)))>0                 
                break;
            end
        end
        t(ind(1),ii) = t(ind(2),jj);
        ind2(i) = ind(2);
    end
end
t(ind2,:) = [];

[p,t]=fixmesh(p,t);




