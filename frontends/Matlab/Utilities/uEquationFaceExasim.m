function [Ru, Rh, B, D, F, G, K, H] = uEquationFaceExasim(dataout, pde, npe, npf, nfe, ne)

% npe = size(mesh.xpe,1);
% npf = size(mesh.xpf,1);
% face = getelemface(pde.nd,pde.elemtype);
% [~,nfe] = size(face);

% ne = size(mesh.t, 2);
%nc = pde.nc;
ncu = pde.ncu;
ncq = pde.ncq;
%nd = pde.nd;

filename = dataout + "uEquationElemFace_B.bin";
tmp = readbin(filename);
B = reshape(tmp, npe, npe, ne, ncu, ncq);

filename = dataout + "uEquationElemFace_D.bin";
tmp = readbin(filename);
D = reshape(tmp, npe, npe, ne, ncu, ncu);

filename = dataout + "uEquationElemFace_Ru.bin";
tmp = readbin(filename);
Ru = reshape(tmp, npe, ne, ncu);

filename = dataout + "uEquationElemFace_F.bin";
tmp = readbin(filename);
F = reshape(tmp, npe, npf*nfe, ne, ncu, ncu);

filename = dataout + "uEquationElemFace_G.bin";
tmp = readbin(filename);
G = reshape(tmp, npf*nfe, npe, ne, ncu, ncq);

filename = dataout + "uEquationElemFace_K.bin";
tmp = readbin(filename);
K = reshape(tmp, npf*nfe, npe, ne, ncu, ncu);

filename = dataout + "uEquationElemFace_H.bin";
tmp = readbin(filename);
H = reshape(tmp, npf*nfe, npf*nfe, ne, ncu, ncu);

filename = dataout + "uEquationElemFace_Rh.bin";
tmp = readbin(filename);
Rh = reshape(tmp, npf*nfe, ne, ncu);





