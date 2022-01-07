function [xg, xx, jac, ug, fg, rne, rqe] = getRqElem(fileout, nge, npe, ncx, nd, ncu, ncq)

fileID = fopen([fileout 'xge.bin'],'r');
xg = fread(fileID,'double'); 
fclose(fileID);
xg = reshape(xg,nge,ncx,[]);

fileID = fopen([fileout 'xxe.bin'],'r');
xx = fread(fileID,'double'); 
fclose(fileID);
xx = reshape(xx,nge,nd*nd,[]);

fileID = fopen([fileout 'jace.bin'],'r');
jac = fread(fileID,'double'); 
fclose(fileID);
jac = reshape(jac,nge,[]);

fileID = fopen([fileout 'uge.bin'],'r');
ug = fread(fileID,'double'); 
fclose(fileID);
ug = reshape(ug,nge,ncu,[]);

fileID = fopen([fileout 'fge.bin'],'r');
fg = fread(fileID,'double'); 
fclose(fileID);
fg = reshape(fg,nge,ncq,[]);

fileID = fopen([fileout 'rne.bin'],'r');
rne = fread(fileID,'double'); 
fclose(fileID);
rne = reshape(rne,npe,ncq,[]);

fileID = fopen([fileout 'rqe.bin'],'r');
rqe = fread(fileID,'double'); 
fclose(fileID);
rqe = reshape(rqe,npe,ncq,[]);

