function mesh1 = extractmesh(mesh, ind, bndexpr)

[p, t] = fixmesh(mesh.p, mesh.t(ind,:));
mesh1 = mkmesh(p, t, mesh.porder, bndexpr, mesh.elemtype, mesh.nodetype);
mesh1.dgnodes = mesh.dgnodes(:,:,ind);



