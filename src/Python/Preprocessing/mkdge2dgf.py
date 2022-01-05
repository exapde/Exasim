from numpy import *

def mkdge2dgf(facecon,entmax):

    #[~,nf] = size(facecon);
    facecon = facecon.flatten('F');
    ndgf = size(facecon);

    ent = sort(unique(facecon));
    ent = ent[ent >-1];
    ndof = size(ent);
    #entmax = ent(end); # largest entity

    # entity-to-index mapping
    ent2ind = -ones(entmax).astype(int);
    ent2ind[ent] = arange(0,ndof);
     # store number of nfighboring elements for each entity
    rowdge2dgf = zeros(ndof).astype(int);
    for i in range(0,ndgf):  # for each face entity
        k = facecon[i];     # entities on element i
        if k>-1:
            ind = ent2ind[k];   # get entity indices
            rowdge2dgf[ind] = rowdge2dgf[ind] + 1;
    rowdge2dgf = concatenate([[0], cumsum(rowdge2dgf)]);

     # store nfighboring-element indices for each entity
    coldge2dgf = zeros(rowdge2dgf[-1]).astype(int);
    inc = zeros(ndof).astype(int);
    for i in range(0,ndgf):  # for each face entity
        k = facecon[i];     # entities on element i
        if k>-1:
            ind = ent2ind[k];   # get entity indices
            # rowdge2dgf(ind): pointer to the list of entities k
            coldge2dgf[rowdge2dgf[ind]+inc[ind]] = i;
            inc[ind] = inc[ind] + 1;

    return rowdge2dgf,coldge2dgf,ent2ind,ent
