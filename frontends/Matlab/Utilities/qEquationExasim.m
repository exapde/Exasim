function [MinvC, MinvE, Mass, Minv, udg, uhat] = qEquationExasim(dataout, pde, npe, npf, nd, ne)

% npe = size(mesh.xpe,1);
% nd = size(mesh.xpe,2);
% npf = size(mesh.xpf,1);
% ne = size(mesh.t, 2);

face = getelemface(nd,pde.elemtype);
[~,nfe] = size(face);

filename = dataout + "qEquationElem_Mass.bin";
tmp = readbin(filename);
Mass = reshape(tmp, npe, npe, ne);

filename = dataout + "qEquationElem_Minv.bin";
tmp = readbin(filename);
Minv = reshape(tmp, npe, npe, ne);

filename = dataout + "qEquationElem_C.bin";
tmp = readbin(filename);
MinvC = reshape(tmp, npe, npe, ne, nd);

filename = dataout + "qEquationFace_E.bin";
tmp = readbin(filename);
MinvE = reshape(tmp, npe, npf*nfe, ne, nd);

filename = dataout + "qEquation_udg.bin";
tmp = readbin(filename);
nc = length(tmp(:))/(npe*ne);
udg = reshape(tmp, npe, nc, ne);

filename = dataout + "qEquation_uhat.bin";
tmp = readbin(filename);
ncu = pde.ncu;
nf = length(tmp(:))/(npf*ncu);
uhat = reshape(tmp, ncu, npf*nf);





