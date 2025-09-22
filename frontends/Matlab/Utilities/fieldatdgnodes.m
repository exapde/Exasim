function UDG1 = fieldatdgnodes(mesh, master, UDG, dgnodes1)

[npe, ~, ne]  = size(dgnodes1);
nc = size(UDG, 2);

[elist, xi] = locatexinmesh(mesh, master, dgnodes1, [], 1e-4);

UDG1 = evalfield(mesh, UDG, elist, xi);

UDG1 = permute(reshape(UDG1, [npe ne nc]),[1 3 2]); 

