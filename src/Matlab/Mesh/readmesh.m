function [p,t] = readmesh(filename, mode)

fileID = fopen(filename,'r');
if mode==0 % binary    
    tmp = fread(fileID,'double');
else % text    
    tmp = fscanf(fileID, '%f');
end
fclose(fileID);

nd = int64(tmp(1));
np = int64(tmp(2));
nve = int64(tmp(3));
ne = int64(tmp(4));
if mode==0 % binary    
    p = reshape(tmp(5:(4+nd*np)),[nd np]);   
    t = reshape(tmp((5+nd*np):(4+nd*np+nve*ne)),[nve ne]);   
else
    p = reshape(tmp(5:(4+nd*np)),[np nd])';   
    t = reshape(tmp((5+nd*np):(4+nd*np+nve*ne)),[ne nve])';   
end

t = int64(t);    


