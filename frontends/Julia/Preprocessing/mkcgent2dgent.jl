function mkcgent2dgent(dgnodes,tol)

dig = Int(log10(1/tol));
dgnodes = round.(dgnodes,digits=dig);

# CG element-to-node connectivities
cgnodes,cgelcon = mkelconcg(dgnodes);

# CG node-to-element connectivities
rowent2elem,colent2elem,~,~,nent = mkent2elem(cgelcon);

cgent2dgent = colent2elem;
npe = size(dgnodes,1);
dim = size(dgnodes,2);
#nent = maximum(cgelcon[:]);

# display(size(dgnodes))
# display(size(cgnodes))
# display(nent)
# display(size(rowent2elem))
for i = 1:nent # for CG node i
    nelem = rowent2elem[i+1]-rowent2elem[i]; # number of elements connected to CG node i
    elem = colent2elem[(rowent2elem[i]+1):rowent2elem[i+1]]; # elements connected to CG node i
    xcg = reshape(cgnodes[:,i],1,dim); # coordinates of CG node i
    for j = 1:nelem # for element j connected to CG node i
        xdg = dgnodes[:,:,elem[j]]; # coordinates of DG nodes on element j
        in = xiny(xcg,xdg,1); # match CG node i to one of the DG nodes on element j
        cgent2dgent[rowent2elem[i]+j] = (elem[j]-1)*npe+in[1]; # index of the matched DG node
    end
end

return cgelcon,rowent2elem,colent2elem,cgent2dgent,cgnodes

end
