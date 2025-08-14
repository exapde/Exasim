function [f,t2fl] = setboundary(p,t,elemtype,bndexpr)

dim = size(p,1);
ne = size(t,2);
face = getelemface(dim,elemtype);
[nvf,nfe] = size(face);
t2fl = reshape(t(face,:),[nvf nfe*ne]);
pf = reshape(p(:,t2fl),[dim nvf nfe ne]);

% interior faces are zero
f = zeros(nfe,ne);

% face-to-element connectivity
f2t = mkf2e(t,elemtype,dim);

% find elements on the domain boundary
ind = find(f2t(3,:)==0);
ind2 = min(2,nvf);

for i = 1:length(ind) % for each element on the domain boundary    
    e = f2t(1,ind(i)); % element e
    l = f2t(2,ind(i)); % local face index
    for k = 1:length(bndexpr) % for each boundary expression
        % evaluate the boundary expression and check if element e belong to this boundary                
        if bndexpr{k}(pf(:,1,l,e)) && bndexpr{k}(pf(:,ind2,l,e)) && bndexpr{k}(pf(:,nvf,l,e)) 
            f(l,e) = k; % then set f(l,e) to k    
            break;
        end
    end
end
