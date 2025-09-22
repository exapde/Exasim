function [tprd, f] = setperiodic(p,t,f,t2fl,prdexpr)

% dim = size(p,1);
% ne = size(t,2);
% face = getelemface(dim,elemtype);
% [nvf,nfe] = size(face);
% t2fl = reshape(t(face,:),[nvf nfe*ne]);

[nfe, ne] = size(f);

f = f(:);
tprd = t;
nprd = size(prdexpr,1);
% periodic boundary faces are negative
for i = 1:nprd % for each periodic expression   
    i1 = f==prdexpr{i,1}; % find faces on the first periodic boundary                 
    f(i1) = -f(i1); % periodic boundary faces are negative
    v1 = t2fl(:,i1); 
    v1 = unique(v1(:)); % vertices on the first periodic boundary
    p1 = prdexpr{i,2}(p(:,v1)); % nodes on the first periodic boundary

    i2 = f==prdexpr{i,3}; % find faces on the second periodic boundary 
    f(i2) = -f(i2); % periodic boundary faces are negative       
    v2 = t2fl(:,i2);
    v2 = unique(v2(:)); % vertices on the second periodic boundary
    p2 = prdexpr{i,4}(p(:,v2)); % nodes on the second periodic boundary
        
    in = xiny(p1',p2'); % match nodes
    v2 = v2(in);        % match vertices

    for j=1:length(v1)
        tprd(tprd==v2(j)) = v1(j);
    end
end

f = reshape(f,[nfe ne]);




