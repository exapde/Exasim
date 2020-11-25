function [f, tprd, t2t, fprd] = facenumbering(p,t,elemtype,bndexpr,prdexpr)

dim = size(p,1);
ne = size(t,2);
face = getelemface(dim,elemtype);
[nvf,nfe] = size(face);

t2fl = reshape(t(face,:),[nvf nfe*ne]);
pf = reshape(p(:,t2fl),[dim nvf nfe ne]);

% interior faces are zero
f = zeros(nfe,ne);

% face-to-element connectivity
[f2t, t2t] = mkf2e(t,elemtype,dim);
% find elements on the domain boundary
ind = find(f2t(3,:)==0);

for i = 1:length(ind) % for each element on the domain boundary    
    e = f2t(1,ind(i)); % element e
    l = f2t(2,ind(i)); % local face index
    for k = 1:length(bndexpr) % for each boundary expression        
        if bndexpr{k}(pf(:,1,l,e)) && bndexpr{k}(pf(:,2,l,e)) && bndexpr{k}(pf(:,nvf,l,e)) % evaluate the boundary expression and check if element e belong to this boundary            
            f(l,e) = k; % then set f(l,e) to k    
            break;
        end
    end
end

pf = reshape(sum(pf,2)/nvf,[dim nfe ne]);
nprd = size(prdexpr,1);
if nprd>0
    f = f(:);
    pf = reshape(pf,[dim nfe*ne]);    
    tprd = t;
    % periodic boundary faces are negative
    fprd = cell(nprd,1);    
    for i = 1:nprd % for each periodic expression   
        i1 = find(f==prdexpr{i,1}); % find faces on the first periodic boundary 
        f(i1) = -f(i1); % periodic boundary faces are negative
        p1 = prdexpr{i,2}(pf(:,i1)); % the centers of faces on the first peiordic boundary        
        e1 = ceil(i1/nfe);      % elements on the first peiordic boundary
        l1 = i1 - (e1-1)*nfe;   % local faces on the first peiordic boundary
        
        i2 = find(f==prdexpr{i,3}); % find faces on the second periodic boundary 
        f(i2) = -f(i2); % periodic boundary faces are negative   
        p2 = prdexpr{i,4}(pf(:,i2)); % the centers of faces on the second peiordic boundary
        e2 = ceil(i2/nfe);      % elements on the second peiordic boundary
        l2 = i2 - (e2-1)*nfe;   % local faces on the second peiordic boundary                
        
        % update t2t to connect periodic elements
        for j=1:length(e1)
            t2t(l2(j),e2(j)) = e1(j);
            t2t(l1(j),e1(j)) = e2(j);
        end
        
        in = xiny(p1',p2');
        fprd{i} = [i1; i2(in)];

        v1 = t2fl(:,i1); 
        v1 = unique(v1(:)); % vertices on the first periodic boundary
        p1 = prdexpr{i,2}(p(:,v1)); % nodes on the first periodic boundary
        
        v2 = t2fl(:,i2(in));
        v2 = unique(v2(:)); % vertices on the second periodic boundary
        p2 = prdexpr{i,4}(p(:,v2)); % nodes on the second periodic boundary
        in = xiny(p1',p2'); % match nodes
        v2 = v2(in);        % match vertices
        
        for j=1:length(v1)
            tprd(tprd==v2(j)) = v1(j);
        end
    end
else
    fprd = [];
    tprd = t;    
end

f = reshape(f,[nfe ne]);




