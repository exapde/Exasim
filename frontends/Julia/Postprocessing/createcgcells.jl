function  createcgcells(cgelcon,tlocal,ne)

nce, nve = size(tlocal);
cells = zeros(Int,ne*nce,nve);
for el=1:ne
    ind = (nce*(el-1)+1):nce*el;
    cells[ind,:] = reshape(cgelcon[tlocal[:],el],(nce,nve));
end
cells = cells .- 1;

return cells

end
