function [AE, FE, DUDG, DUDG_DUH] = uEquationSchurExasim2(dataout, pde, npe, npf, nfe, ne)

ncu = pde.ncu;

filename = dataout + "uEquationElemSchur_Ru.bin";
tmp = readbin(filename);
DUDG = reshape(tmp, npe, ncu, ne);

filename = dataout + "uEquationElemSchur_Rh.bin";
tmp = readbin(filename);
FE = reshape(tmp, ncu, npf*nfe, ne);

filename = dataout + "uEquationElemSchur_DinvF.bin";
tmp = readbin(filename);
DUDG_DUH = reshape(tmp, npe, ncu, ncu, npf*nfe, ne);

filename = dataout + "uEquationElemSchur_DinvH.bin";
tmp = readbin(filename);
AE = reshape(tmp, ncu, npf*nfe, ncu, npf*nfe, ne);



