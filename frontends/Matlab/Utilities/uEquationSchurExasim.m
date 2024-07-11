function [AE, FE, DUDG, DUDG_DUH, D, F, K, H] = uEquationSchurExasim(dataout, pde, npe, npf, nfe, ne)

% npe = size(mesh.xpe,1);
% npf = size(mesh.xpf,1);
% face = getelemface(pde.nd,pde.elemtype);
% [~,nfe] = size(face);

%ne = size(mesh.t, 2);
ncu = pde.ncu;

filename = dataout + "uEquationElemSchur_D.bin";
tmp = readbin(filename);
D = reshape(tmp, npe, ncu, npe, ncu, length(tmp)/(npe*ncu*npe*ncu));

filename = dataout + "uEquationElemSchur_F.bin";
tmp = readbin(filename);
F = reshape(tmp, npe, ncu, ncu, npf*nfe, length(tmp)/(npe*ncu*npf*nfe*ncu));

filename = dataout + "uEquationElemSchur_K.bin";
tmp = readbin(filename);
K = reshape(tmp, ncu, npf*nfe, npe, ncu, length(tmp)/(npe*ncu*npf*nfe*ncu));

filename = dataout + "uEquationElemSchur_H.bin";
tmp = readbin(filename);
H = reshape(tmp, ncu, npf*nfe, ncu, npf*nfe, length(tmp)/(npf*nfe*ncu*npf*nfe*ncu));

filename = dataout + "uEquationElemSchur_Ru.bin";
tmp = readbin(filename);
DUDG = reshape(tmp, npe, ncu, length(tmp)/(npe*ncu));

filename = dataout + "uEquationElemSchur_Rh.bin";
tmp = readbin(filename);
FE = reshape(tmp, ncu, npf*nfe, length(tmp)/(ncu*npf*nfe));

filename = dataout + "uEquationElemSchur_DinvF.bin";
tmp = readbin(filename);
DUDG_DUH = reshape(tmp, npe, ncu, ncu, npf*nfe, length(tmp)/(npe*ncu*npf*nfe*ncu));

filename = dataout + "uEquationElemSchur_DinvH.bin";
tmp = readbin(filename);
AE = reshape(tmp, ncu, npf*nfe, ncu, npf*nfe, length(tmp)/(npf*nfe*ncu*npf*nfe*ncu));



