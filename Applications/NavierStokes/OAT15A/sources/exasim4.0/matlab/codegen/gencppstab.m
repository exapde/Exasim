function gencppstab(stab,xdg,uhg,udg1,udg2,odg1,odg2,nlg,param,time,appname)

% void FhatDriver(dstype *fg, dstype *xg, dstype *uh, dstype *ug1, dstype *ug2, dstype * og1, 
%      dstype * og2, dstype *nl, meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, 
%      tempstruct &tmp, commonstruct &common, Int ngf, Int f1, Int f2, Int backend)

filename = ['opuStab' appname '.cpp'];
%delete(filename);
gid = fopen(filename,'w');

for d = 1:length(stab)
if ~isempty(stab{d})    
    
ncu = numel(stab{d});    
ccode(stab{d}(:),'file','tmp.c');

fid  = fopen('tmp.c','r');
f=fread(fid,'*char')';
fclose(fid);
f = strrep(f, 't0', 'A0[0][0]');    
fid  = fopen('tmp.c','w');
fprintf(fid,'%s',f);
fclose(fid);

str = ['template <typename T> void opuStab' appname num2str(d) 'd(T *f, T *xdg, T *uhg, T *udg1, T *udg2, T *odg1, T *odg2, T *nlg, T *param, T time, int ng, int nc, int ncu, int nd, int ncx, int nco)'];
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
for i=1:length(uhg{d})
    str = ['T uhg' num2str(i) ' = uhg[' num2str(i-1) '*ng+i];'];
    fprintf(gid, '\t\t%s\n', str);                  
end
for i=1:length(udg1{d})
    str = ['T udg1' num2str(i) ' = udg1[' num2str(i-1) '*ng+i];'];
    fprintf(gid, '\t\t%s\n', str);                  
end
for i=1:length(udg2{d})
    str = ['T udg2' num2str(i) ' = udg2[' num2str(i-1) '*ng+i];'];
    fprintf(gid, '\t\t%s\n', str);                  
end
for i=1:length(odg1{d})
    str = ['T odg1' num2str(i) ' = odg1[' num2str(i-1) '*ng+i];'];
    fprintf(gid, '\t\t%s\n', str);                  
end
for i=1:length(odg2{d})
    str = ['T odg2' num2str(i) ' = odg2[' num2str(i-1) '*ng+i];'];
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
            %strj = ['f[' num2str(j) '*ng+i] = 0.0;'];
            %fprintf(gid, '\t\t%s\n', strj);                  
        end
        a1 = a2+1;              
    end

    str = strrep(str, '  ', '');
    str = strrep(str, 'A0[', 'f[');
    str = strrep(str, '][0] =', '*ng+i] +=');                          
    if isempty(i1)==1
        str = ['double ' str];
    end

    fprintf(gid, '\t\t%s\n', str);    
    tline = fgetl(fid);        
    i=i+1;   
end
% if a1<ncu
%     for j = a1:(ncu-1)                
%         strj = ['f[' num2str(j) '*ng+i] = 0.0;'];
%         fprintf(gid, '\t\t%s\n', strj);                  
%     end
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

str = ['template <typename T> void opuStab' appname '(T *f, T *xdg, T *uhg, T *udg1, T *udg2, T *odg1, T *odg2, T *nlg, T *param, T time, int ng, int nc, int ncu, int nd, int ncx, int nco)'];
fprintf(gid, '%s\n', str); 
str = '{';
fprintf(gid, '%s\n', str);         
str = 'if (nd == 1) {';
fprintf(gid, '\t%s\n', str);         
if ~isempty(stab{1})
str = ['opuStab' appname num2str(1) 'd(f, xdg, uhg, udg1, udg2, odg1, odg2, nlg, param, time, ng, nc, ncu, nd, ncx, nco);'];
fprintf(gid, '\t\t%s\n', str);      
end
str = '}';
fprintf(gid, '\t%s\n', str);         
str = 'else if (nd == 2) {';
fprintf(gid, '\t%s\n', str);         
if ~isempty(stab{2})
str = ['opuStab' appname num2str(2) 'd(f, xdg, uhg, udg1, udg2, odg1, odg2, nlg, param, time, ng, nc, ncu, nd, ncx, nco);'];
fprintf(gid, '\t\t%s\n', str);      
end
str = '}';
fprintf(gid, '\t%s\n', str);         
str = 'else {';
fprintf(gid, '\t%s\n', str);         
if ~isempty(stab{3})
str = ['opuStab' appname num2str(3) 'd(f, xdg, uhg, udg1, udg2, odg1, odg2, nlg, param, time, ng, nc, ncu, nd, ncx, nco);'];
fprintf(gid, '\t\t%s\n', str);      
end
str = '}';
fprintf(gid, '\t%s\n', str);         
str = '}';
fprintf(gid, '%s\n', str);             
fprintf(gid, '\n'); 

str = ['template void opuStab' appname '(double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int);'];
fprintf(gid, '%s\n', str);  
str = ['template void opuStab' appname '(float *, float *, float *, float *, float *, float *, float *, float *, float *, float,  int, int, int, int, int, int);'];
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

filename2 = ['cpuStab' appname '.cpp'];
fid  = fopen(filename2,'w');
fprintf(fid,'%s',f);
fclose(fid);

% CUDA code
fid  = fopen(filename,'r');
f=fread(fid,'*char')';
fclose(fid);

f = strrep(f, 'opu', 'gpu');    
f = strrep(f, 'template <typename T> void', 'template <typename T>  __global__  void');    
f = strrep(f, ['__global__  void gpuStab' appname '(T'], ['void gpuStab' appname '(T']);   

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

oldstr = '(f, xdg, uhg,';
str = '<<<gridDim, blockDim>>>(f, xdg, uhg,';
f = strrep(f, oldstr, str);

ind1 = strfind(f,['void gpuStab' appname '(T']);
ind2 = strfind(f,'i += blockDim.x * gridDim.x;');
for i = length(ind2):-1:1
    if ind2(i)>ind1        
        f(ind2(i):ind2(i)+30) = '';
    end
end

filename2 = ['gpuStab' appname '.cu'];
fid  = fopen(filename2,'w');
fprintf(fid,'%s',f);
fclose(fid);



