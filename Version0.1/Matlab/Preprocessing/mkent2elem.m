function [rowent2elem,colent2elem,ent2ind,ent] = mkent2elem(elcon)

[~,ne] = size(elcon);
ent = unique(elcon(:));
ent = ent(ent>0);
ndof = length(ent);
entmax = ent(end); % largest entity

% entity-to-index mapping
ent2ind = zeros(entmax,1);
ent2ind(ent) = (1:ndof);

 % store number of neighboring elements for each entity
rowent2elem = zeros(ndof,1);
for i = 1:ne  % for each element i
    elc = elcon(:,i);     % entities on element i  
    k   = unique(elc(:)); % remove duplicate entities on \element i  
    k   = k(k>0);
    ind = ent2ind(k);   % get entity indices 
    rowent2elem(ind) = rowent2elem(ind) + 1;
end
rowent2elem = [0; cumsum(rowent2elem)];

 % store neighboring-element indices for each entity
colent2elem = zeros(rowent2elem(end),1); 
inc = ones(ndof,1);
for i = 1:ne
    elc = elcon(:,i);   % entities on element i  
    k = unique(elc(:)); % remove duplicate entities on element i  
    k   = k(k>0);
    ind = ent2ind(k);   % get entity indices 
    % rowent2elem(ind): pointer to the list of entities k
    colent2elem(rowent2elem(ind)+inc(ind)) = i;
    inc(ind) = inc(ind) + 1;
end

