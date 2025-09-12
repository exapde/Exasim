function [v_new, C_new] = removenodes(v, C, j_list)
% Removes multiple nodes v(j_list) from the node list v and updates connectivity matrix C
% v      : 1 x N list of node values
% C      : M x N connectivity matrix
% j_list : indices (into v) of nodes to delete

    N = numel(v);
    j_list = unique(j_list(:)', 'stable');  % ensure row vector
    if any(j_list > N) || any(j_list < 1)
        error("j_list contains invalid indices");
    end

    % Node values to remove
    v_removed = zeros(1,N);

    % Step 1: Map node value to column index
    max_node = max(v);
    idx_of = zeros(1, max_node);
    for i = 1:N
        idx_of(v(i)) = i;
    end

    % Step 2: Initialize
    v_new = v;
    C_new = C;
    
    % Step 3: Update all neighbor columns of deleted nodes
    for j = j_list        
        vj = v(j);  % node to be removed
        v_removed(j) = vj;
        neighbors = C_new(2 : 1 + C_new(1,j), j);
        num_neighbors = numel(neighbors);        
        
        for k = 1:num_neighbors
            vi = neighbors(k);  % neighbor of vj
            if ~ismember(vi, v_removed)  % only update if vi is not being deleted
                i = idx_of(vi);

                % Remove vj from vi's neighbor list
                ni = C_new(1,i);
                nbrs = C_new(2:1+ni, i);
                nbrs(nbrs == vj) = [];

                % Add other neighbors of vj (excluding vi and deleted nodes)
                new_links = setdiff(neighbors, [vi, v_removed]);
                nbrs = union(nbrs, new_links);  % avoid duplicates

                % Update column i
                C_new(:,i) = 0;
                C_new(1,i) = numel(nbrs);
                C_new(2:1+numel(nbrs), i) = nbrs;
            end
        end        
    end

    % Step 4: Remove columns of deleted nodes
    keep_mask = true(1, N);
    keep_mask(j_list) = false;
    v_new = v_new(keep_mask);
    C_new = C_new(:, keep_mask);
end
