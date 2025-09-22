function R = matvecExasim(pde, mesh)

npe = size(mesh.xpe,1);
npf = size(mesh.xpf,1);
face = getelemface(pde.nd,pde.elemtype);
[~,nfe] = size(face);

ne = size(mesh.t, 2);
nf = size(mesh.f2t,2);
ncu = pde.ncu;

dataout = pde.buildpath + "/dataout/";

filename = dataout + "outhdgMatVec.bin";
tmp = readbin(filename);
R = reshape(tmp, ncu, npf*nf);



