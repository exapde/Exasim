function mkdge2dgf(facecon,entmax)

#[~,nf] = size(facecon);
ndgf = length(facecon);

ent = sort(unique(facecon[:]));
ent = ent[ent .>0 ];
ndof = length(ent);
#entmax = ent(end); # largest entity

# entity-to-index mapping
ent2ind = zeros(Int,entmax,1);
ent2ind[ent] = 1:ndof;
 # store number of nfighboring elements for each entity
rowdge2dgf = zeros(Int,ndof,1);
for i = 1:ndgf  # for each face entity
    k = facecon[i];     # entities on element i
    if k>0
        ind = ent2ind[k];   # get entity indices
        rowdge2dgf[ind] = rowdge2dgf[ind] + 1;
    end
end
rowdge2dgf=[0; cumsum(rowdge2dgf[:])];

 # store nfighboring-element indices for each entity
coldge2dgf = zeros(Int,rowdge2dgf[end],1);
inc = ones(Int,ndof,1);
for i = 1:ndgf
    k = facecon[i];     # entities on element i
    if k>0
        ind = ent2ind[k];   # get entity indices
        # rowdge2dgf(ind): pointer to the list of entities k
        coldge2dgf[rowdge2dgf[ind]+inc[ind]] = i;
        inc[ind] = inc[ind] + 1;
    end
end

return rowdge2dgf,coldge2dgf,ent2ind,ent

end
