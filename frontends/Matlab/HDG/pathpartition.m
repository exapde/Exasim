function [paths,nelems,npaths] = pathpartition(e2e, K)

    N = size(e2e,2);
    visited = false(1, N);
    paths = zeros(1,N);
    nelems = zeros(1,N+1);
    path = zeros(1,K);
    npaths = 0;    
    nelem = 0;
    nep = 0;
    last_node = -1;  % Track last node of the previous path    
    
    while nelem < N
        if last_node == -1
            % No prior path: pick the first unvisited node
            start_node = find(~visited, 1);            
        else
            % Try to extend from last_node
            neighbors_raw = e2e(:, last_node);
            neighbors = neighbors_raw(neighbors_raw > 0);
            candidate_starts = neighbors(~visited(neighbors));
            
            kk = 1;
            while isempty(candidate_starts) && nep>kk     
              last_node = path(end-kk);
              neighbors_raw = e2e(:, last_node);
              neighbors = neighbors_raw(neighbors_raw > 0);
              candidate_starts = neighbors(~visited(neighbors));              
              kk = kk + 1;
            end
              
            if isempty(candidate_starts)
                % No neighbor to extend from, fallback to first unvisited
                start_node = find(~visited, 1);
            else
                % Pick candidate with fewest unvisited neighbors
                scores1 = zeros(1, length(candidate_starts));
                scores2 = zeros(1, length(candidate_starts));
                for i = 1:length(candidate_starts)
                    ni = candidate_starts(i);
                    ni_neighbors_raw = e2e(:, ni);
                    ni_neighbors = ni_neighbors_raw(ni_neighbors_raw > 0);
                    scores1(i) = sum(~visited(ni_neighbors));
                    scores2(i) = length(ni_neighbors);
                end
                %[~, idx] = min(scores);
                idx = find(scores1==min(scores1));
                if length(idx)>1
                  [~, i2] = min(scores2(idx));
                  idx = idx(i2);                  
                end                
                
                start_node = candidate_starts(idx);
            end
        end
        
        path = grow_flexible_path(start_node, visited, e2e, K);
        nep = length(path);
        if nep<K
          path = grow_flexible_path(path(end), visited, e2e, K);
        end        
        visited(path) = true;        
        paths(nelem+1:nelem+nep) = path;
        last_node = path(end);                
        npaths = npaths + 1;
        nelem = nelem + nep;
        nelems(1+npaths) = nelem;
    end
    nelems = nelems(1:npaths+1);
end

function path = grow_flexible_path(start_node, visited, e2e, K)

    path = zeros(1,K);
    path(1) = start_node;
    current = start_node;
    kk = 1;
    while true
        neighbors_raw = e2e(:, current);            % Get all M entries
        neighbors = neighbors_raw(neighbors_raw > 0); % Remove padding zeros
        unvisited_neighbors = neighbors(~visited(neighbors));
        unvisited_neighbors = setdiff(unvisited_neighbors, path); % avoid cycles

        if isempty(unvisited_neighbors)
            break;
        end

        % Greedy: choose neighbor with fewest unvisited neighbors

        % Compute score for each unvisited neighbor (number of unvisited neighbors it has)
        scores1 = zeros(1, length(unvisited_neighbors));
        scores2 = zeros(1, length(unvisited_neighbors));
        for i = 1:length(unvisited_neighbors)
            ni = unvisited_neighbors(i);
            ni_neighbors_raw = e2e(:, ni);
            ni_neighbors = ni_neighbors_raw(ni_neighbors_raw > 0); % Remove padding zeros
            scores1(i) = sum(~visited(ni_neighbors));
            scores2(i) = length(ni_neighbors);
        end
        
        % Pick the neighbor with the fewest unvisited neighbors (greedy)
        idx = find(scores1==min(scores1));
        if length(idx)>1
          [~, i2] = min(scores2(idx));
          idx = idx(i2);                  
        end                
        next_node = unvisited_neighbors(idx);

        kk = kk + 1;
        path(kk) = next_node;
        current = next_node;

        if kk >= K
            break;
        end
    end
    path = path(1:kk);
end
