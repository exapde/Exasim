function tmp = readbin(filename)

fileID = fopen(filename,'r');
tmp = fread(fileID,'double');
fclose(fileID);

