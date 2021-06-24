function mesh = updatemesh(mesh,bcm)

% --------------- update MESH structure ----------------%
% Reorder faces
%[mesh.f,mesh.t2f,mesh.f2f,mesh.flev] = faceordering(mesh.f, mesh.t2f); 
[facecon,f,t2f] = mkconnectivity(mesh.p, mesh.t, mesh.porder, mesh.elemtype, mesh.bndexpr, []);
%figure(2);clf;plotface(mesh.p, mesh.t, f);

[mesh.f,mesh.t2f,mesh.facecon] = reorderface(f,t2f,facecon,bcm,mesh.perm);
%figure(3);clf;plotface(mesh.p, mesh.t, mesh.f);

mesh.bf = reshape(mesh.f(abs(mesh.t2f'),end),[size(mesh.perm,2) size(mesh.t,1)]);
mesh.f2f = mkf2f(mesh.f, mesh.t2f);

