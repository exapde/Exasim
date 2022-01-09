function createhighordermesh(mesh,app)

display("run createnodes...");
#[mesh.dgnodes,mesh.elemtype,perm] = createnodes(mesh.p,mesh.t,app.porder);
mesh.dgnodes, elemtype, perm = createnodes(mesh.p,mesh.t,app.porder);

display("run facenumbering...");
mesh.f, mesh.tprd = facenumbering(mesh.p,mesh.t,elemtype,mesh.boundaryexpr,mesh.periodicexpr);

# Project nodes on the curved boundaries
if length(mesh.curvedboundaryexpr)>0 && app.porder>1 && maximum(mesh.curvedboundary[:])>0
    display("Project dgnodes onto the curved boundaries...");
    fd = mesh.curvedboundaryexpr;
    nd = size(mesh.p,1);
    nfe,ne = size(mesh.f);
    if nd == 2
        for i = 1:ne
            for j = 1:nfe
                if mesh.f[j,i]!=0 # boundary element
                    k = abs(mesh.f[j,i]); # get boundary index
                    if mesh.curvedboundary[k]==1 # if this boundary is curved
                        p = mesh.dgnodes[perm[:,j],:,i];
                        deps = 1e-8*maximum(maximum(p,dim=1)-minimum(p,dim=1));
                        d = fd[k](p);
                        dgradx = (fd[k]([p[:,1]+deps p[:,2]])-d)/deps;
                        dgrady = (fd[k]([p[:,1] p[:,2]+deps])-d)/deps;
                        dgrad2 = dgradx.*dgradx + dgrady.*dgrady;
                        dgrad2[dgrad2==0] = 1;
                        p = p-[d.*dgradx./dgrad2 d.*dgrady./dgrad2];
                        mesh.dgnodes[perm[:,j],:,i] = p;
                    end
                end
            end
        end
    elseif nd==3
        for i = 1:ne
            for j = 1:nfe
                if mesh.f[j,i]!=0 # boundary element
                    k = abs(mesh.f[j,i]); # get boundary index
                    if mesh.curvedboundary[k]==1 # if this boundary is curved
                        p = mesh.dgnodes[perm[:,j],:,i];
                        deps = 1e-8*maximum(maximum(p,dim=1)-minimum(p,dim=1));
                        d = fd[k](p);
                        dgradx = (fd[k]([p[:,1]+deps,p[:,2],p[:,3]])-d)/deps;
                        dgrady = (fd[k]([p[:,1],p[:,2]+deps,p[:,3]])-d)/deps;
                        dgradz = (fd[k]([p[:,1],p[:,2],p[:,3]+deps])-d)/deps;
                        dgrad2 = dgradx.*dgradx+dgrady.*dgrady+dgradz.*dgradz;
                        dgrad2[dgrad2==0] = 1;
                        p = p-[d.*dgradx./dgrad2 d.*dgrady./dgrad2 d.*dgradz./dgrad2];
                        mesh.dgnodes[perm[:,j],:,i] = p;
                    end
                end
            end
        end
    end
end

return mesh

end
