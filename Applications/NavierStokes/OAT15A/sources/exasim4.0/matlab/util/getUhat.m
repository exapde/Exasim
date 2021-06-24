function [xg, nlg, jac, ug1, ug2, uh] = getUhat(fileout, npf, ncx, nd, ncu, nc)

fileID = fopen([fileout 'xgh.bin'],'r');
xg = fread(fileID,'double'); 
fclose(fileID);
xg = reshape(xg,npf,ncx,[]);

fileID = fopen([fileout 'nlgh.bin'],'r');
nlg = fread(fileID,'double'); 
fclose(fileID);
nlg = reshape(nlg,npf,nd,[]);

fileID = fopen([fileout 'jach.bin'],'r');
jac = fread(fileID,'double'); 
fclose(fileID);
jac = reshape(jac,npf,[]);

fileID = fopen([fileout 'ug1h.bin'],'r');
ug1 = fread(fileID,'double'); 
fclose(fileID);
ug1 = reshape(ug1,npf,nc,[]);

fileID = fopen([fileout 'ug2h.bin'],'r');
ug2 = fread(fileID,'double'); 
fclose(fileID);
ug2 = reshape(ug2,npf,nc,[]);

fileID = fopen([fileout 'uh.bin'],'r');
uh = fread(fileID,'double'); 
fclose(fileID);
uh = reshape(uh,npf,ncu,[]);




