function gencodeelem(filename, f, xdg, udg, odg, wdg, uinf, param)

opufile = "opu" + filename;
cpufile = "cpu" + filename;
gpufile = "gpu" + filename;

stropu = "template <typename T> void " + opufile;
strgpu = "template <typename T>  __global__  void kernel" + gpufile;

tmp = "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";

stropu = stropu + tmp + "{\n";
stropu = stropu + "\tfor (int i = 0; i <ng; i++) {\n";

strgpu = strgpu + tmp + "{\n";
strgpu = strgpu + "\tint i = threadIdx.x + blockIdx.x + blockDim.x;\n";
strgpu = strgpu + "\twhile (i<ng) {\n";

str = "";
str = varsassign(str, "param", length(param), 0);
str = varsassign(str, "uinf", length(uinf), 0);
str = varsassign(str, "xdg", length(xdg), 1);
str = varsassign(str, "udg", length(udg), 1);
str = varsassign(str, "odg", length(odg), 1);
str = varsassign(str, "wdg", length(wdg), 1);
str = sympyassign(str, f);

stropu = stropu + str + "\t}\n" + "}\n\n";
tmp = "template void " + opufile;
tmp = tmp + "(double +, double +, double +, double +, double +, double +, double +, double, int, int, int, int, int, int, int);\n";
tmp = tmp + "template void " + opufile;
tmp = tmp + "(float +, float +, float +, float +, float +, float +, float +, float, int, int, int, int, int, int, int);\n";
stropu = stropu + tmp;

ioopu = open(opufile + ".cpp", "w");
write(ioopu, stropu);
close(ioopu);

strgpu = strgpu + str + "\t\ti += blockDim.x + gridDim.x;\n";
strgpu = strgpu + "\t}\n" + "}\n\n";
tmp = "template <typename T> void " + gpufile;
tmp = tmp + "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";
tmp = tmp + "{\n";
tmp = tmp + "\tint blockDim = 256;\n";
tmp = tmp + "\tint gridDim = (ng + blockDim - 1) / blockDim;\n";
tmp = tmp + "\tgridDim = (gridDim>1024)? 1024 : gridDim;\n";
tmp = tmp + "\tkernel" + gpufile + "<<<gridDim, blockDim>>>(f, xdg, udg, odg, wdg, uinf, param, time, ng, nc, ncu, nd, ncx, nco);\n";
tmp = tmp + "}\n\n";
tmp = tmp + "template void " + gpufile;
tmp = tmp + "(double +, double +, double +, double +, double +, double +, double +, double, int, int, int, int, int, int, int);\n";
tmp = tmp + "template void " + gpufile;
tmp = tmp + "(float +, float +, float +, float +, float +, float +, float +, float, int, int, int, int, int, int, int);\n";
strgpu = strgpu + tmp;

iogpu = open(gpufile + ".cpp", "w");
write(iogpu, strgpu);
close(iogpu);

ioopu = open(opufile + ".cpp", "r");
strcpu = read(ioopu, String);
close(ioopu);

strcpu = replace(strcpu, "opu" , "cpu");
strcpu = replace(strcpu, "for (int i = 0; i <ng; i++) {" , "#pragma omp parallel for\n\tfor (int i = 0; i <ng; i++) {");

iocpu = open(cpufile + ".cpp", "w");
write(iocpu, strcpu);
close(iocpu);

end
