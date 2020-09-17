

[MiC, MiE] = massinvnd(master, mesh.dgnodes);

fileID = fopen('app/MiC2.bin','r');
MiCt = fread(fileID,'double');
fclose(fileID);
MiCt = reshape(MiCt,size(MiC));
max(max(abs(MiC(:)+MiCt(:))))

fileID = fopen('app/MiE.bin','r');
MiEt = fread(fileID,'double');
fclose(fileID);
MiEt = reshape(MiEt,size(MiE));
max(max(abs(MiE(:)-MiEt(:))))
