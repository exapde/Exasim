function [Ru, udg, uh] = getResidual(fileout, npe, npf, ncu, nc)

fileID = fopen([fileout '_uh.bin'],'r');
uh = fread(fileID,'double'); 
fclose(fileID);
uh = reshape(uh,npf,ncu,[]);

fileID = fopen([fileout '_udg.bin'],'r');
udg = fread(fileID,'double'); 
fclose(fileID);
udg = reshape(udg,npe,nc,[]);

fileID = fopen([fileout '_Ru.bin'],'r');
Ru = fread(fileID,'double'); 
fclose(fileID);
Ru = reshape(Ru,npe,ncu,[]);
