from numpy import *
from getelemface import *

def  mkv2t(t,ne):

    ndof = max(t.flatten())+1;
    re = zeros(ndof).astype(int); # store number of neighboring elements for each entity
    for i in range(0,ne):
        k = t[:,i];
        k = k[k>-1];
        re[k] = re[k] + 1;

    re = hstack((array([0]), cumsum(re)));
    re = re.astype(int);

    ce = zeros(re[-1]).astype(int);  # store neighboring-element indices for each entity
    ind = ones(ndof).astype(int);
    for i in range(0,ne):
        k = t[:,i];   # entities on element i
        k = k[k>-1];
        # re(k): pointer to the list of entities k
        ce[re[k]+ind[k]-1] = i+1;
        ind[k] = ind[k] + 1;

    return re,ce
