function gencodeelem2(filename::String, f, xdg, udg, odg, wdg, uinf, param, time)

foldername = "app";
opufile = "opu" * filename;
cpufile = "cpu" * filename;
gpufile = "gpu" * filename;

stropu = "template <typename T> void " * opufile;
strgpu = "template <typename T>  __global__  void kernel" * gpufile;

tmp = "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne)\n";

stropu = stropu * tmp * "{\n";
stropu = stropu * "\tfor (int i = 0; i <ng; i++) {\n";

strgpu = strgpu * tmp * "{\n";
strgpu = strgpu * "\tint i = threadIdx.x + blockIdx.x * blockDim.x;\n";
strgpu = strgpu * "\twhile (i<ng) {\n";

str = "\t\tint j = i%npe;\n";
str = str * "\t\tint k = (i-j)/npe;\n";

fstr = string(f[:]);       
str = varsassign(str, "param", length(param), 0, fstr);
str = varsassign(str, "uinf", length(uinf), 0, fstr);

varname = "xdg";
for i = 1:length(xdg)
    str1 = varname * string(i);
    str2 = varname * "[j+npe*" * string(i-1) * "+npe*ncx*k" * "]";
    str = str * "\t\tT " * str1 * " = " * str2 * ";\n";
end
varname = "udg";
for i = 1:length(udg)
    str1 = varname * string(i);
    str2 = varname * "[j+npe*" * string(i-1) * "+npe*nc*k" * "]";
    str = str * "\t\tT " * str1 * " = " * str2 * ";\n";
end
varname = "odg";
for i = 1:length(odg)
    str1 = varname * string(i);
    str2 = varname * "[j+npe*" * string(i-1) * "+npe*nco*k" * "]";
    str = str * "\t\tT " * str1 * " = " * str2 * ";\n";
end
varname = "wdg";
for i = 1:length(wdg)
    str1 = varname * string(i);
    str2 = varname * "[j+npe*" * string(i-1) * "+npe*ncw*k" * "]";
    str = str * "\t\tT " * str1 * " = " * str2 * ";\n";
end

n = length(f);
#for i=1:n
#    f[i] = sympy.simplify(f[i]);
#end
ts, fs = sympy.cse(f);
m = length(ts);
for i = 1:m
    str1 = sympy.ccode(ts[i][1]);
    str2 = sympy.ccode(ts[i][2]);
    str = str * "\t\tT " * str1 * " = " * str2 * ";\n";
end
for i = 1:n
    str1 = "\t\tf[j+npe*" * string(i-1) * "+npe*nce*k" * "]";
    str2 = sympy.ccode(fs[1][i]);
    str = str * str1 * " = " * str2 * ";\n";
end

stropu = stropu * str * "\t}\n" * "}\n\n";
tmp = "template void " * opufile;
tmp = tmp * "(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int, int, int);\n";
tmp = tmp * "template void " * opufile;
tmp = tmp * "(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int, int, int);\n";
stropu = stropu * tmp;

ioopu = open(foldername * "/" * opufile * ".cpp", "w");
write(ioopu, stropu);
close(ioopu);

strgpu = strgpu * str * "\t\ti += blockDim.x * gridDim.x;\n";
strgpu = strgpu * "\t}\n" * "}\n\n";
tmp = "template <typename T> void " * gpufile;
tmp = tmp * "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne)\n";
tmp = tmp * "{\n";
tmp = tmp * "\tint blockDim = 256;\n";
tmp = tmp * "\tint gridDim = (ng + blockDim - 1) / blockDim;\n";
tmp = tmp * "\tgridDim = (gridDim>1024)? 1024 : gridDim;\n";
tmp = tmp * "\tkernel" * gpufile * "<<<gridDim, blockDim>>>(f, xdg, udg, odg, wdg, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne);\n";
tmp = tmp * "}\n\n";
tmp = tmp * "template void " * gpufile;
tmp = tmp * "(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int, int, int);\n";
tmp = tmp * "template void " * gpufile;
tmp = tmp * "(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int, int, int);\n";
strgpu = strgpu * tmp;

iogpu = open(foldername * "/" * gpufile * ".cu", "w");
write(iogpu, strgpu);
close(iogpu);

strcpu = replace(stropu, "opu" => "cpu");
strcpu = replace(strcpu, "for (int i = 0; i <ng; i++) {" => "#pragma omp parallel for\n\tfor (int i = 0; i <ng; i++) {");

iocpu = open(foldername * "/" * cpufile * ".cpp", "w");
write(iocpu, strcpu);
close(iocpu);

end
