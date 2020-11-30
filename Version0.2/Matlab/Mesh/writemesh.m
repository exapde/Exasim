function writemesh(p, t, filename, mode)

endian = 'native';
[nd,np] = size(p);
[nve,ne] = size(t);
tmp = [nd np nve ne];

fileID = fopen(filename,'w');
if mode==0 % binary    
    fwrite(fileID,tmp,'double',endian);
    fwrite(fileID,p,'double',endian);
    fwrite(fileID,t,'double',endian);    
else % text    
    fprintf(fileID,'%d %d %d %d\n',tmp);
    if nd==1
        fprintf(fileID,'%.16f\n',p);
    elseif nd==2
        fprintf(fileID,'%.16f %.16f\n',p);
    elseif nd==3
        fprintf(fileID,'%.16f %.16f %.16f\n',p);
    end
    if nve==2
        fprintf(fileID,'%d %d\n',t);
    elseif nve==3
        fprintf(fileID,'%d %d %d\n',t);        
    elseif nve==4        
        fprintf(fileID,'%d %d %d %d\n',t);        
    elseif nve==8
        fprintf(fileID,'%d %d %d %d %d %d %d %d\n',t);        
    end    
end
fclose(fileID);

