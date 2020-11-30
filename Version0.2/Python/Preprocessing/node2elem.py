from numpy import *

def node2elem(nodes,re,ce):

    # # node-2-elemsent connectivity
    # [re,ce] = mkv2t(t,ne);

    # elemsent-to-node connectivity in nonoverlapping subdomain Omega_i
    # nodes = t(:,intelem);

    # remove duplications
    nodes = sort(unique(nodes.flatten(order='F')));

    # remove zero or negative
    nodes = nodes[nodes>-1];

    # find all elemsents connected to nodes
    nt = sum(re[nodes+1]-re[nodes]);
    elems = zeros(nt).astype(int);
    j = 0;
    for i in range(0,nodes.size):
        ni = nodes[i];
        ri = arange(re[ni]+1,re[ni+1]+1);
        elems[j:(j+ri.size)] = ce[ri-1];
        j = j+ri.size;

    elems = sort(unique(elems))-1;

    return elems
