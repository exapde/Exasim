using Preprocessing

function createcggrid(dgnodes,tlocal)

# Get some dimensions of the mesh
~, nd, ne = size(dgnodes);
~, nve    = size(tlocal);

cgelcon,~,~,~,cgnodes = Preprocessing.mkcgent2dgent(dgnodes,1e-7);

cgnodes = cgnodes';
cgcells = createcgcells(cgelcon,tlocal,ne);

elemtype = 0;
if (nd==2) && (nve==4)
    elemtype=1;
end
if (nd==3) && (nve==8)
    elemtype=1;
end

cell_t = getcelltype(nd,elemtype);

return cgnodes, cgelcon, cgcells, cell_t

end
