function [a1,a2] = comparebinfiles(fn1, fn2)

fileID = fopen([fn1 '.bin'],'r');
a1 = fread(fileID,'double');
fclose(fileID);

fileID = fopen([fn2 '.bin'],'r');
a2 = fread(fileID,'double');
fclose(fileID);

    