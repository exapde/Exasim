function [group_id, group_graph] = group_nodes_with_two_neighbors(ind, f2f)
% C: M-by-N connectivity matrix
% Output:
%   group_id: 1-by-N array assigning each node to a group
%   group_graph: adjacency cell array describing connections between groups

nbf = findneighborfaces(ind, f2f);

p = size(f2f,1);
n = length(ind);
group_id = -ones(1, n); 
group_elem = zeros(p, n); 
group_graph = zeros(2,n);

fij = zeros(2*p,p);
nij = zeros(1,p);

group_counter = 0;
for i = 1:n
  fi = ind(i); % face fi
  j = nbf(1,i); % number of neighbors of face fi  
  if j<=2
    group_counter = group_counter + 1;
    group_id(i) = group_counter;
    group_elem(1,group_counter) = 1;
    group_elem(2,group_counter) = fi;    
%     k = group_graph(1,group_counter) + 1;
%     group_graph(1,group_counter) = k;
  else    
    nbfi = nbf(2:(j+1),i); % neighbors of face fi        
    for m = 1:j     % loop over each neighbor of face fi      
      fj = nbfi(m); % face fj
      l = ind==fj;
      k = nbf(1,l);      
      nbfj = nbf(2:(k+1),l); % neighbors of face fj
      % neighbors of both fi and fj      
      nij(m) = 0;
      for i1 = 1:j
        if i1 ~= m
          nij(m) = nij(m)+1;
          fij(nij(m),m) = nbfi(i1);
        end
      end
      for i1 = 1:k
        if ismember(nbfj(i1), nbfi) == 0
          nij(m) = nij(m)+1;
          fij(nij(m),m) = nbfj(i1);
        end
      end            
    end
    
  end
end

end

% function fij = join_neighbors(nbfi, nbfj, j, k, fi, fj)
% 
% 
% end

% [M, N] = size(C);
% group_id = -ones(1, N);  % initialize group assignment
% group_counter = 0;
% 
% % Build adjacency list
% neighbors = cell(1, N);
% for i = 1:N
%     num_nbrs = C(1, i);
%     neighbors{i} = C(2:(1+num_nbrs), i)';
% end
% 
% group_graph = containers.Map('KeyType','int32', 'ValueType','any');
% 
% % Group assignment via BFS-like traversal
% for i = 1:N
%     if group_id(i) ~= -1
%         continue;
%     end
% 
%     group_counter = group_counter + 1;
%     q = i;
%     group_id(i) = group_counter;
% 
%     while ~isempty(q)
%         u = q(1);
%         q(1) = [];
% 
%         for v = neighbors{u}
%             if group_id(v) == -1
%                 group_id(v) = group_counter;
%                 q(end+1) = v; %#ok<AGROW>
%             elseif group_id(v) ~= group_id(u)
%                 gu = group_id(u);
%                 gv = group_id(v);
%                 % Record inter-group connection
%                 if ~isKey(group_graph, gu)
%                     group_graph(gu) = [];
%                 end
%                 if ~isKey(group_graph, gv)
%                     group_graph(gv) = [];
%                 end
%                 if ~ismember(gv, group_graph(gu))
%                     group_graph(gu) = [group_graph(gu), gv];
%                 end
%                 if ~ismember(gu, group_graph(gv))
%                     group_graph(gv) = [group_graph(gv), gu];
%                 end
%             end
%         end
%     end
% end
% 
% % Warn if any group connects to more than 2 other groups
% keys_g = keys(group_graph);
% for i = 1:length(keys_g)
%     g = keys_g{i};
%     neighbors_g = group_graph(g);
%     if length(neighbors_g) > 2
%         warning('Group %d connects to more than 2 neighboring groups. Refinement may be needed.', g);
%     end
% end
% end
