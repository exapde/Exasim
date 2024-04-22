function BE = blockjacobiExasim(dataout, pde, npf, nf)

% npe = size(mesh.xpe,1);
% npf = size(mesh.xpf,1);
% face = getelemface(pde.nd,pde.elemtype);
% [~,nfe] = size(face);

% ne = size(mesh.t, 2);
%nf = size(mesh.f2t,2);
ncu = pde.ncu;
ncf = ncu*npf;

filename = dataout + "hdgBlockJacobi.bin";
tmp = readbin(filename);
BE = reshape(tmp, ncf, ncf, nf);



