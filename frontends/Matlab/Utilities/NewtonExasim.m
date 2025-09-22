function [AE, FE, DUDG, DUDG_DUH] = NewtonExasim(dataout, pde, npe, npf, nfe)

% npe = size(mesh.xpe,1);
% npf = size(mesh.xpf,1);
% face = getelemface(pde.nd,pde.elemtype);
% [~,nfe] = size(face);

%ne = size(mesh.t, 2);
ncu = pde.ncu;

filename = dataout + "newton_DUDG.bin";
tmp = readbin(filename);
DUDG = reshape(tmp, npe, ncu, length(tmp)/(npe*ncu));

filename = dataout + "newton_FE.bin";
tmp = readbin(filename);
FE = reshape(tmp, ncu, npf*nfe, length(tmp)/(ncu*npf*nfe));

filename = dataout + "newton_DUDG_DUH.bin";
tmp = readbin(filename);
DUDG_DUH = reshape(tmp, npe, ncu, ncu, npf*nfe, length(tmp)/(npe*ncu*npf*nfe*ncu));

filename = dataout + "newton_AE.bin";
tmp = readbin(filename);
AE = reshape(tmp, ncu, npf*nfe, ncu, npf*nfe, length(tmp)/(npf*nfe*ncu*npf*nfe*ncu));



