function stropu = gencodeelem2(filename, f, xdg, udg, odg, wdg, uinf, param, time, foldername)

%foldername = "app";
opufile = "opu" + filename;
cpufile = "cpu" + filename;
gpufile = "gpu" + filename;

stropu = "template <typename T> void " + opufile;
strgpu = "template <typename T>  __device__  void device" + gpufile;

tmp = "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne)\n";

stropu = stropu + tmp + "{\n";
stropu = stropu + "\tfor (int i = 0; i <ng; i++) {\n";

strgpu = strgpu + tmp + "{\n";
strgpu = strgpu + "\tint i = threadIdx.x + blockIdx.x * blockDim.x;\n";
strgpu = strgpu + "\twhile (i<ng) {\n";

str = "\t\tint j = i%%npe;\n";
str = str + "\t\tint k = (i-j)/npe;\n";

fstr = string(f(:));
str = varsassign(str, "param", length(param), 0, fstr);
str = varsassign(str, "uinf", length(uinf), 0, fstr);

varname = "xdg";
for i = 1:length(xdg)
    str1 = varname + string(i);
    str2 = varname + "[j+npe*" + string(i-1) + "+npe*ncx*k" + "]";
    str = str + "\t\tT " + str1 + " = " + str2 + ";\n";
end
varname = "udg";
for i = 1:length(udg)
    str1 = varname + string(i);
    str2 = varname + "[j+npe*" + string(i-1) + "+npe*nc*k" + "]";
    str = str + "\t\tT " + str1 + " = " + str2 + ";\n";
end
varname = "odg";
for i = 1:length(odg)
    str1 = varname + string(i);
    str2 = varname + "[j+npe*" + string(i-1) + "+npe*nco*k" + "]";
    str = str + "\t\tT " + str1 + " = " + str2 + ";\n";
end
varname = "wdg";
for i = 1:length(wdg)
    str1 = varname + string(i);
    str2 = varname + "[j+npe*" + string(i-1) + "+npe*ncw*k" + "]";
    str = str + "\t\tT " + str1 + " = " + str2 + ";\n";
end

n = length(f);
% for i=1:n
%     f(i) = sympy.simplify(f(i));
% end
ccode(f(:),'file','tmp.c');

fid  = fopen('tmp.c','r');
f=fread(fid,'*char')';
fclose(fid);
f = strrep(f, 't0', 'A0[0][0]');    
fid  = fopen('tmp.c','w');
fprintf(fid,'%s',f);
fclose(fid);

mystr = str;
fid = fopen('tmp.c','r'); 
tline = fgetl(fid); 
i=1; a1 = 0;       
while ischar(tline)        
    str = tline;

    i1 = strfind(str,'[');        
    i2 = strfind(str,']');        
    if isempty(i1)==0    
        a2 = str2num(str((i1(1)+1):(i2(1)-1)));                        
        for j = a1:(a2-1)                
            strj = ['f[j+npe*' num2str(j) '+npe*nce*k] = 0.0;'];
            mystr = mystr + "\t\t" + string(strj) + "\n";
        end
        a1 = a2+1;              
    end

    str = strrep(str, '  ', '');
    str = strrep(str, 'A0[', 'f[j+npe*');
    str = strrep(str, '][0]', '+npe*nce*k]');                          
    if isempty(i1)==1
        str = "T " + string(str);
    end

    mystr = mystr + "\t\t" + string(str) + "\n";
    tline = fgetl(fid);        
    i=i+1;   
end
if a1<n
    for j = a1:(n-1)                
        strj = ['f[j+npe*' num2str(j) '+npe*nce*k] = 0.0;'];
        mystr = mystr + "\t\t" + string(strj) + "\n";
    end
end
fclose(fid);
str = mystr;

stropu = stropu + str + "\t}\n" + "}\n\n";
tmp = "template void " + opufile;
tmp = tmp + "(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int, int, int);\n";
tmp = tmp + "template void " + opufile;
tmp = tmp + "(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int, int, int);\n";
stropu = stropu + tmp;

fid = fopen(foldername + "/" + opufile + ".cpp", 'w');
fprintf(fid, char(stropu));
fclose(fid);

strgpu = strgpu + str + "\t\ti += blockDim.x * gridDim.x;\n";
strgpu = strgpu + "\t}\n" + "}\n\n";

tmp = "template <typename T> __global__ void kernel" + gpufile;
tmp = tmp + "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne)\n";
tmp = tmp + "{\n";
tmp = tmp + "\tdevice" + gpufile + "(f, xdg, udg, odg, wdg, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne);\n";
tmp = tmp + "}\n\n";

tmp = tmp + "template <typename T> void " + gpufile;
tmp = tmp + "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne)\n";
tmp = tmp + "{\n";
tmp = tmp + "\tint blockDim = 256;\n";
tmp = tmp + "\tint gridDim = (ng + blockDim - 1) / blockDim;\n";
tmp = tmp + "\tgridDim = (gridDim>1024)? 1024 : gridDim;\n";
tmp = tmp + "\tkernel" + gpufile + "<<<gridDim, blockDim>>>(f, xdg, udg, odg, wdg, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne);\n";
tmp = tmp + "}\n\n";
tmp = tmp + "template void " + gpufile;
tmp = tmp + "(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int, int, int);\n";
tmp = tmp + "template void " + gpufile;
tmp = tmp + "(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int, int, int);\n";
strgpu = strgpu + tmp;

fid = fopen(foldername + "/" + gpufile + ".cu", 'w');
fprintf(fid, char(strgpu));
fclose(fid);

strcpu = strrep(stropu, 'opu', "cpu");
strcpu = strrep(strcpu, "for (int i = 0; i <ng; i++) {", "#pragma omp parallel for\n\tfor (int i = 0; i <ng; i++) {");

iocpu = fopen(foldername + "/" + cpufile + ".cpp", 'w');
fprintf(fid, char(strcpu));
fclose(iocpu);

delete(char("tmp.c"));

end
