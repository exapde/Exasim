from numpy import *

def mkent2elem(elcon):

    ne = elcon.shape[1]
    ent = unique(elcon.flatten('F'));
    ent = ent[ent>-1];
    ndof = ent.size;
    entmax = ent[-1]; # largest entity

    # entity-to-index mapping
    ent2ind = zeros(entmax+1).astype(int);
    ent2ind[ent] = arange(0,ndof);

     # store number of neighboring elements for each entity
    rowent2elem = zeros(ndof).astype(int);
    for i in range(0,ne):  # for each element i
        elc = elcon[:,i];     # entities on element i
        k   = unique(elc.flatten('F')); # remove duplicate entities on \element i
        k   = k[k>-1];
        ind = ent2ind[k];   # get entity indices
        rowent2elem[ind] = rowent2elem[ind] + 1;

    rowent2elem = concatenate([[0], cumsum(rowent2elem.flatten())]);

     # store neighboring-element indices for each entity
    colent2elem = zeros(int(rowent2elem[-1])).astype(int);
    inc = zeros(ndof).astype(int);
    for i in range(0,ne):
        elc = elcon[:,i];   # entities on element i
        k = unique(elc[:]); # remove duplicate entities on element i
        k   = k[k>-1];
        ind = ent2ind[k];   # get entity indices
        # rowent2elem(ind): pointer to the list of entities k
        colent2elem[rowent2elem[ind]+inc[ind]] = i;
        inc[ind] = inc[ind] + 1;
        
    return rowent2elem,colent2elem,ent2ind,ent
