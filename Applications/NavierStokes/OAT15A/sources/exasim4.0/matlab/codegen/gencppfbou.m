function gencppfbou(fbou,xdg,udg,uhg,odg,nlg,tau,uinf,param,time,appname)

filename = ['opuFbou' appname '.cpp'];
%delete(filename);
gid = fopen(filename,'w');

for d = 1:length(fbou)
if ~isempty(fbou{d})
[ncu, nbc] = size(fbou{d});
for k = 1:nbc

ccode(fbou{d}(:,k),'file','tmp.c');
fid  = fopen('tmp.c','r');
f=fread(fid,'*char')';
fclose(fid);
f = strrep(f, 't0', 'A0[0][0]');    
fid  = fopen('tmp.c','w');
fprintf(fid,'%s',f);
fclose(fid);

str = ['template <typename T> void opuFbou' appname num2str(d) 'd' num2str(k) '(T *f, T *xdg, T *udg, T *uhg, T *odg, T *nlg, T *tau, T *uinf, T *param, T time, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco)'];
fprintf(gid, '%s\n', str); 
str = '{';
fprintf(gid, '%s\n', str);         

str = 'for (int i = 0; i <ng; i++) {';
fprintf(gid, '\t%s\n', str);         

for i=1:length(param{d})
    str = ['T param' num2str(i) ' = param[' num2str(i-1) '];'];
    fprintf(gid, '\t\t%s\n', str);                  
end        
for i=1:length(uinf{d})
    str = ['T uinf' num2str(i) ' = uinf[' num2str(i-1) '];'];
    fprintf(gid, '\t\t%s\n', str);                  
end        
for i=1:length(tau{d})
    str = ['T tau' num2str(i) ' = tau[' num2str(i-1) '];'];
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
for i=1:length(uhg{d})
    str = ['T uhg' num2str(i) ' = uhg[' num2str(i-1) '*ng+i];'];
    fprintf(gid, '\t\t%s\n', str);                  
end
for i=1:length(odg{d})
    str = ['T odg' num2str(i) ' = odg[' num2str(i-1) '*ng+i];'];
    fprintf(gid, '\t\t%s\n', str);                  
end
for i=1:length(nlg{d})
    str = ['T nlg' num2str(i) ' = nlg[' num2str(i-1) '*ng+i];'];
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

fclose(fid);

str = '}';
fprintf(gid, '\t%s\n', str);             

str = '}';
fprintf(gid, '%s\n', str);             
fprintf(gid, '\n'); 

delete('tmp.c');

end
end
end

str = ['template <typename T> void opuFbou' appname '(T *f, T *xdg, T *udg, T *uhg, T *odg, T *nlg, T *tau, T *uinf, T *param, T time, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco)'];
fprintf(gid, '%s\n', str); 
str = '{';
fprintf(gid, '%s\n', str);         
str = 'if (nd == 1) {';
fprintf(gid, '\t%s\n', str);         
if ~isempty(fbou{1})
[~, nbc] = size(fbou{1});  
for k = 1:nbc
    if k == 1
        str = ['if (ib == ' num2str(k) ')'];    
    elseif k<nbc
        str = ['else if (ib == ' num2str(k) ')'];    
    else            
        str = ['else if (ib == ' num2str(k) ')'];    
    end 
    fprintf(gid, '\t\t%s\n', str);      
    str = ['opuFbou' appname num2str(1) 'd' num2str(k) '(f, xdg, udg, uhg, odg, nlg, tau, uinf, param, time, ib, ng, nc, ncu, nd, ncx, nco);'];
    fprintf(gid, '\t\t\t%s\n', str);                  
end
end
str = '}';
fprintf(gid, '\t%s\n', str);         
str = 'else if (nd == 2) {';
fprintf(gid, '\t%s\n', str);    
if ~isempty(fbou{2})
[~, nbc] = size(fbou{2});  
for k = 1:nbc
    if k == 1
        str = ['if (ib == ' num2str(k) ')'];    
    elseif k<nbc
        str = ['else if (ib == ' num2str(k) ')'];    
    else            
        str = ['else if (ib == ' num2str(k) ')'];    
    end 
    fprintf(gid, '\t\t%s\n', str);      
    str = ['opuFbou' appname num2str(2) 'd' num2str(k) '(f, xdg, udg, uhg, odg, nlg, tau, uinf, param, time, ib, ng, nc, ncu, nd, ncx, nco);'];
    fprintf(gid, '\t\t\t%s\n', str);                  
end
end
str = '}';
fprintf(gid, '\t%s\n', str);         
str = 'else {';
fprintf(gid, '\t%s\n', str);     
if ~isempty(fbou{3})
[~, nbc] = size(fbou{3});  
for k = 1:nbc
    if k == 1
        str = ['if (ib == ' num2str(k) ')'];    
    elseif k<nbc
        str = ['else if (ib == ' num2str(k) ')'];    
    else            
        str = ['else if (ib == ' num2str(k) ')'];    
    end 
    fprintf(gid, '\t\t%s\n', str);      
    str = ['opuFbou' appname num2str(3) 'd' num2str(k) '(f, xdg, udg, uhg, odg, nlg, tau, uinf, param, time, ib, ng, nc, ncu, nd, ncx, nco);'];
    fprintf(gid, '\t\t\t%s\n', str);                  
end
end
str = '}';
fprintf(gid, '\t%s\n', str);         
str = '}';
fprintf(gid, '%s\n', str);             
fprintf(gid, '\n'); 

str = ['template void opuFbou' appname '(double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int);'];
fprintf(gid, '%s\n', str);  
str = ['template void opuFbou' appname '(float *, float *, float *, float *, float *, float *, float *, float *, float *, float,  int, int, int, int, int, int, int);'];
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

filename2 = ['cpuFbou' appname '.cpp'];
fid  = fopen(filename2,'w');
fprintf(fid,'%s',f);
fclose(fid);

% CUDA code
fid  = fopen(filename,'r');
f=fread(fid,'*char')';
fclose(fid);

f = strrep(f, 'opu', 'gpu');    
f = strrep(f, 'template <typename T> void', 'template <typename T>  __global__  void');    
f = strrep(f, ['__global__  void gpuFbou' appname '(T'], ['void gpuFbou' appname '(T']);   

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

ind1 = strfind(f,['void gpuFbou' appname '(T']);
ind2 = strfind(f,'i += blockDim.x * gridDim.x;');
for i = length(ind2):-1:1
    if ind2(i)>ind1        
        f(ind2(i):ind2(i)+30) = '';
    end
end

filename2 = ['gpuFbou' appname '.cu'];
fid  = fopen(filename2,'w');
fprintf(fid,'%s',f);
fclose(fid);



