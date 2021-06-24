function gencpptdfunc(tdfunc,xdg,udg,odg,param,time,appname)

filename = ['opuTdfunc' appname '.cpp'];
%delete(filename);
gid = fopen(filename,'w');

for d = 1:length(tdfunc)
if ~isempty(tdfunc{d})    
    
ncu = numel(tdfunc{d});    
ccode(tdfunc{d}(:),'file','tmp.c');

fid  = fopen('tmp.c','r');
f=fread(fid,'*char')';
fclose(fid);
f = strrep(f, 't0', 'A0[0][0]');    
fid  = fopen('tmp.c','w');
fprintf(fid,'%s',f);
fclose(fid);

str = ['template <typename T> void opuTdfunc' appname num2str(d) 'd(T *f, T *xdg, T *udg, T *odg, T *param, T time, int ng, int nc, int ncu, int nd, int ncx, int nco)'];
fprintf(gid, '%s\n', str); 
str = '{';
fprintf(gid, '%s\n', str);         

str = 'for (int i = 0; i <ng; i++) {';
fprintf(gid, '\t%s\n', str);         

for i=1:length(param{d})
    str = ['T param' num2str(i) ' = param[' num2str(i-1) '];'];
    fprintf(gid, '\t\t%s\n', str);                  
end        

for i=1:length(xdg{d})
    str = ['T xdg' num2str(i) ' = xdg[' num2str(i-1) '*ng+i];'];
    fprintf(gid, '\t\t%s\n', str);                  
end
for i=1:length(udg{d})
    str = ['T udg' num2str(i) ' = udg[' num2str(i-1) '*ng+i];'];
    fprintf(gid, '\t\t%s\n', str);                  
end
for i=1:length(odg{d})
    str = ['T odg' num2str(i) ' = odg[' num2str(i-1) '*ng+i];'];
    fprintf(gid, '\t\t%s\n', str);                  
end
fprintf(gid, '\n');      

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
            strj = ['f[' num2str(j) '*ng+i] = 0.0;'];
            fprintf(gid, '\t\t%s\n', strj);                  
        end
        a1 = a2+1;              
    end

    str = strrep(str, '  ', '');
    str = strrep(str, 'A0[', 'f[');
    str = strrep(str, '][0]', '*ng+i]');                          
    if isempty(i1)==1
        str = ['double ' str];
    end

    fprintf(gid, '\t\t%s\n', str);    
    tline = fgetl(fid);        
    i=i+1;   
end
if a1<ncu
    for j = a1:(ncu-1)                
        strj = ['f[' num2str(j) '*ng+i] = 0.0;'];
        fprintf(gid, '\t\t%s\n', strj);                  
    end
end

% tline = fgetl(fid); 
% n = 0; 
% while ischar(tline)        
%     str = tline;        
%     i0 = strfind(str,'A0');        
%     if (isempty(i0)==1)
%         str = ['f[' num2str(n) '*ng+i] = 0.0;'];
%         fprintf(gid, '\t\t%s\n', str);    
%     else
%         str = strrep(str, '  ', '');
%         str = strrep(str, 'A0[', 'f[');
%         str = strrep(str, '][0]', '*ng+i]');       
%         fprintf(gid, '\t\t%s\n', str);    
%     end
%     tline = fgetl(fid);        
%     n = n + 1; 
% end
% for j = n:(ncu*d-1)                
%     str = ['f[' num2str(n) '*ng+i] = 0.0;'];
%     fprintf(gid, '\t\t%s\n', str);                  
% end

fclose(fid);

str = '}';
fprintf(gid, '\t%s\n', str);             

str = '}';
fprintf(gid, '%s\n', str);             
fprintf(gid, '\n'); 

delete('tmp.c');

end
end

str = ['template <typename T> void opuTdfunc' appname '(T *f, T *xdg, T *udg, T *odg, T *param, T time, int ng, int nc, int ncu, int nd, int ncx, int nco)'];
fprintf(gid, '%s\n', str); 
str = '{';
fprintf(gid, '%s\n', str);         
str = 'if (nd == 1) {';
fprintf(gid, '\t%s\n', str);     
if ~isempty(tdfunc{1})
str = ['opuTdfunc' appname num2str(1) 'd(f, xdg, udg, odg, param, time, ng, nc, ncu, nd, ncx, nco);'];
fprintf(gid, '\t\t%s\n', str);      
end
str = '}';
fprintf(gid, '\t%s\n', str);         
str = 'else if (nd == 2) {';
fprintf(gid, '\t%s\n', str);         
if ~isempty(tdfunc{2})
str = ['opuTdfunc' appname num2str(2) 'd(f, xdg, udg, odg, param, time, ng, nc, ncu, nd, ncx, nco);'];
fprintf(gid, '\t\t%s\n', str);      
end
str = '}';
fprintf(gid, '\t%s\n', str);         
str = 'else {';
fprintf(gid, '\t%s\n', str);         
if ~isempty(tdfunc{3})
str = ['opuTdfunc' appname num2str(3) 'd(f, xdg, udg, odg, param, time, ng, nc, ncu, nd, ncx, nco);'];
fprintf(gid, '\t\t%s\n', str);      
end
str = '}';
fprintf(gid, '\t%s\n', str);         
str = '}';
fprintf(gid, '%s\n', str);             
fprintf(gid, '\n'); 

str = ['template void opuTdfunc' appname '(double *, double *, double *, double *, double *, double, int, int, int, int, int, int);'];
fprintf(gid, '%s\n', str);  
str = ['template void opuTdfunc' appname '(float *, float *, float *, float *, float *, float,  int, int, int, int, int, int);'];
fprintf(gid, '%s\n', str);  
fprintf(gid, '\n'); 

fclose(gid);            

% OpenMP code
fid  = fopen(filename,'r');
f=fread(fid,'*char')';
fclose(fid);

f = strrep(f, 'opu', 'cpu');    
str = ['#pragma omp parallel for' char(13) char(10) char(9) 'for (int i = 0; i <ng; i++) {'];
f = strrep(f, 'for (int i = 0; i <ng; i++) {', str);    

filename2 = ['cpuTdfunc' appname '.cpp'];
fid  = fopen(filename2,'w');
fprintf(fid,'%s',f);
fclose(fid);

% CUDA code
fid  = fopen(filename,'r');
f=fread(fid,'*char')';
fclose(fid);

f = strrep(f, 'opu', 'gpu');    
f = strrep(f, 'template <typename T> void', 'template <typename T>  __global__  void');    
f = strrep(f, ['__global__  void gpuTdfunc' appname '(T'], ['void gpuTdfunc' appname '(T']);   

str = ['int i = threadIdx.x + blockIdx.x * blockDim.x;' char(13) char(10) char(9) 'while (i<ng) {'];
f = strrep(f, 'for (int i = 0; i <ng; i++) {', str);    

oldstr = [char(9) '}'];
str = [char(10) char(9) char(9) 'i += blockDim.x * gridDim.x;' char(13) char(10) char(9) '}'];
f = strrep(f, oldstr, str);    

oldstr = 'if (nd == 1)';
str = ['int blockDim = 256;' char(10) char(9) ...
       'int gridDim = (ng + blockDim - 1) / blockDim;' char(10) char(9) ...
       'gridDim = (gridDim>1024)? 1024 : gridDim;' char(10) char(9)...
       'if (nd == 1)'];
f = strrep(f, oldstr, str);

oldstr = '(f, xdg, udg,';
str = '<<<gridDim, blockDim>>>(f, xdg, udg,';
f = strrep(f, oldstr, str);

ind1 = strfind(f,['void gpuTdfunc' appname '(T']);
ind2 = strfind(f,'i += blockDim.x * gridDim.x;');
for i = length(ind2):-1:1
    if ind2(i)>ind1        
        f(ind2(i):ind2(i)+30) = '';
    end
end

filename2 = ['gpuTdfunc' appname '.cu'];
fid  = fopen(filename2,'w');
fprintf(fid,'%s',f);
fclose(fid);



