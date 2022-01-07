function [Mass, Minv] = getMassInv(fileout)

fileID = fopen([fileout 'Mass.bin'],'r');
Mass = fread(fileID,'double'); 
fclose(fileID);

fileID = fopen([fileout 'Minv.bin'],'r');
Minv = fread(fileID,'double'); 
fclose(fileID);

