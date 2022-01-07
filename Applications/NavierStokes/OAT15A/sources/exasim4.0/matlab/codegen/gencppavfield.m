function gencppavfield(avfield,xdg,udg,odg,param,time,appname)

filename = ['opuAVfield' appname '.cpp'];
%delete(filename);
gid = fopen(filename,'w');

for d = 1:length(avfield)
if ~isempty(avfield{d})    
    
ccode(avfield{d}(:),'file','tmp.c');

fid  = fopen('tmp.c','r');
f=fread(fid,'*char')';
fclose(fid);
f = strrep(f, 't0', 'A0[0][0]');    
fid  = fopen('tmp.c','w');
fprintf(fid,'%s',f);
fclose(fid);

str = ['template <typename T> void opuAVfield' appname num2str(d) 'd(T *f, T *xdg, T *udg, T *odg, T *param, T time, int ng, int nc, int ncu, int nd, int ncx, int nco, int npe, int ne)'];
fprintf(gid, '%s\n', str); 
str = '{';
fprintf(gid, '%s\n', str);         

str = 'for (int i = 0; i <ng; i++) {';
fprintf(gid, '\t%s\n', str);         

str = 'int j = i%npe;';
fprintf(gid, '\t\t%s\n', str);         
str = 'int k = (i-j)/npe;';
fprintf(gid, '\t\t%s\n', str);         

for i=1:length(param{d})
    str = ['T param' num2str(i) ' = param[' num2str(i-1) '];'];
    fprintf(gid, '\t\t%s\n', str);                  
end        

for i=1:length(xdg{d})
    str = ['T xdg' num2str(i) ' = xdg[j+npe*' num2str(i-1) '+npe*ncx*k];'];
    fprintf(gid, '\t\t%s\n', str);                  
end
for i=1:length(udg{d})
    str = ['T udg' num2str(i) ' = udg[j+npe*' num2str(i-1) '+npe*nc*k];'];
    fprintf(gid, '\t\t%s\n', str);                  
end
for i=1:length(odg{d})
    str = ['T odg' num2str(i) ' = odg[j+npe*' num2str(i-1) '+npe*nco*k];'];
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
            strj = ['f[j+npe*' num2str(j) '+npe*nco*k] = 0.0;'];
            fprintf(gid, '\t\t%s\n', strj);                  
        end
        a1 = a2+1;              
    end

    str = strrep(str, '  ', '');
    str = strrep(str, 'A0[', 'f[j+npe*');
    str = strrep(str, '][0]', '+npe*nco*k]');                          
    if isempty(i1)==1
        str = ['double ' str];
    end

    fprintf(gid, '\t\t%s\n', str);    
    tline = fgetl(fid);        
    i=i+1;   
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

str = ['template <typename T> void opuAVfield' appname '(T *f, T *xdg, T *udg, T *odg, T *param, T time, int ng, int nc, int ncu, int nd, int ncx, int nco, int npe, int ne)'];
fprintf(gid, '%s\n', str); 
str = '{';
fprintf(gid, '%s\n', str);         
str = 'if (nd == 1) {';
fprintf(gid, '\t%s\n', str);         
if ~isempty(avfield{1})
str = ['opuAVfield' appname num2str(1) 'd(f, xdg, udg, odg, param, time, ng, nc, ncu, nd, ncx, nco, npe, ne);'];
fprintf(gid, '\t\t%s\n', str);      
end
str = '}';
fprintf(gid, '\t%s\n', str);         
str = 'else if (nd == 2) {';
fprintf(gid, '\t%s\n', str);         
if ~isempty(avfield{2})
str = ['opuAVfield' appname num2str(2) 'd(f, xdg, udg, odg, param, time, ng, nc, ncu, nd, ncx, nco, npe, ne);'];
fprintf(gid, '\t\t%s\n', str);      
end
str = '}';
fprintf(gid, '\t%s\n', str);         
str = 'else {';
fprintf(gid, '\t%s\n', str);         
if ~isempty(avfield{3})
str = ['opuAVfield' appname num2str(3) 'd(f, xdg, udg, odg, param, time, ng, nc, ncu, nd, ncx, nco, npe, ne);'];
fprintf(gid, '\t\t%s\n', str);      
end
str = '}';
fprintf(gid, '\t%s\n', str);         
str = '}';
fprintf(gid, '%s\n', str);             
fprintf(gid, '\n'); 

str = ['template void opuAVfield' appname '(double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int);'];
fprintf(gid, '%s\n', str);  
str = ['template void opuAVfield' appname '(float *, float *, float *, float *, float *, float,  int, int, int, int, int, int, int, int);'];
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

filename2 = ['cpuAVfield' appname '.cpp'];
fid  = fopen(filename2,'w');
fprintf(fid,'%s',f);
fclose(fid);

% CUDA code
fid  = fopen(filename,'r');
f=fread(fid,'*char')';
fclose(fid);

f = strrep(f, 'opu', 'gpu');    
f = strrep(f, 'template <typename T> void', 'template <typename T>  __global__  void');    
f = strrep(f, ['__global__  void gpuAVfield' appname '(T'], ['void gpuAVfield' appname '(T']);   

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

ind1 = strfind(f,['void gpuAVfield' appname '(T']);
ind2 = strfind(f,'i += blockDim.x * gridDim.x;');
for i = length(ind2):-1:1
    if ind2(i)>ind1        
        f(ind2(i):ind2(i)+30) = '';
    end
end

filename2 = ['gpuAVfield' appname '.cu'];
fid  = fopen(filename2,'w');
fprintf(fid,'%s',f);
fclose(fid);



