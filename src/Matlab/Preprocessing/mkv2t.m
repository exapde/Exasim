function [re,ce] = mkv2t(t,ne)

sz = size(t);
if sz(1)==ne
    t = t';
end

ndof = max(t(:));
re = zeros(ndof,1); % store number of neighboring elements for each entity
for i = 1:ne
    k = t(:,i);
    k = k(k>0);
    re(k) = re(k) + 1;
end
re=[0; cumsum(re)];

ce = zeros(re(end),1);  % store neighboring-element indices for each entity
in = ones(ndof,1);
for i = 1:ne
    k = t(:,i);   % entities on element i  
    k = k(k>0);
    % re(k): pointer to the list of entities k
    ce(re(k)+in(k)) = i;
    in(k) = in(k) + 1;
end
