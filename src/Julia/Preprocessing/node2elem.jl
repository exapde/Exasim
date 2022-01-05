function node2elem(nodes,re,ce)

# # node-2-elemsent connectivity
# [re,ce] = mkv2t(t,ne);

# elemsent-to-node connectivity in nonoverlapping subdomain Omega_i
# nodes = t(:,intelem);

# remove duplications
nodes = sort(unique(nodes[:]));

# remove zero or negative
nodes = nodes[nodes.>0];

# find all elemsents connected to nodes
nt = sum(re[nodes.+1]-re[nodes]);
elems = zeros(Int,nt,1);
j = 1;
for i = 1:length(nodes)
    ni = nodes[i];
    ri = (re[ni]+1):re[ni+1];
    elems[j:(j+length(ri)-1)] = ce[ri];
    j = j+length(ri);
end
elems = sort(unique(elems));

return elems

end
