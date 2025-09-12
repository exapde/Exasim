function tmp = readbin(filename, type)

if nargin < 2
  fileID = fopen(filename,'r');
  tmp = fread(fileID,'double');
  fclose(fileID);
else
  fileID = fopen(filename,'r');
  tmp = fread(fileID,type);
  fclose(fileID);
end

