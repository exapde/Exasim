function mystr = getccode(f, varstr)

mystr = "";

n = length(f);
ccode(f(:),'file','tmp.c');

fid  = fopen('tmp.c','r');
f=fread(fid,'*char')';
fclose(fid);
f = strrep(f, 't0', 'A0[0][0]');    
fid  = fopen('tmp.c','w');
fprintf(fid,'%s',f);
fclose(fid);

fid = fopen('tmp.c','r'); 
tline = fgetl(fid); 
i=1; a1 = 0;       
while ischar(tline)        
    str = tline;

    i1 = strfind(str,'[');        
    i2 = strfind(str,']');        
    if isempty(i1)==0    
        a2 = str2num(str((i1+1):(i2-1)));                        
        for j = a1:(a2-1)                
            strj = [varstr num2str(j) '*ng+i] = 0.0;'];
            mystr = mystr + "\t\t" + string(strj) + "\n";
        end
        a1 = a2+1;              
    end

    str = strrep(str, '  ', '');
    str = strrep(str, 'A0[', varstr);
    str = strrep(str, '][0]', '*ng+i]');                          
    if isempty(i1)==1
        str = "T " + string(str);
    end

    mystr = mystr + "\t\t" + string(str) + "\n";
    tline = fgetl(fid);        
    i=i+1;   
end
if a1<n
    for j = a1:(n-1)                
        strj = [varstr num2str(j) '*ng+i] = 0.0;'];
        mystr = mystr + "\t\t" + string(strj) + "\n";
    end
end
fclose(fid);

delete(char("tmp.c"));

end
