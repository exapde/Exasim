function writebin(filename,a)

endian = 'native';
fileID = fopen(filename,'w');
fwrite(fileID,a(:),'double',endian);
fclose(fileID);



