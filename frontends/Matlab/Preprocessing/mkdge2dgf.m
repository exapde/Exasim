function [rowdge2dgf,coldge2dgf,ent2ind,ent] = mkdge2dgf(facecon,entmax)

%[~,nf] = size(facecon);
ndgf = numel(facecon);

ent = unique(facecon(:));
ent = ent(ent>0);
ndof = length(ent);
%entmax = ent(end); % largest entity

% entity-to-index mapping
ent2ind = zeros(entmax,1);
ent2ind(ent) = (1:ndof);

 % store number of neighboring elements for each entity
rowdge2dgf = zeros(ndof,1);
for i = 1:ndgf  % for each face entity
    k = facecon(i);     % entities on element i
    if k>0
        ind = ent2ind(k);   % get entity indices 
        rowdge2dgf(ind) = rowdge2dgf(ind) + 1;
    end
end

rowdge2dgf=[0; cumsum(rowdge2dgf)];
 % store neighboring-element indices for each entity
coldge2dgf = zeros(rowdge2dgf(end),1); 
inc = ones(ndof,1);
for i = 1:ndgf % for every face entity
    k = facecon(i);      
    if k>0
        ind = ent2ind(k);   % get entity indices 
        % rowdge2dgf(ind): pointer to the list of entities k
        coldge2dgf(rowdge2dgf(ind)+inc(ind)) = i;
        inc(ind) = inc(ind) + 1;
    end
end

