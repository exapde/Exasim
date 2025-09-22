function stropu = gencodeface3(filename, f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, foldername)

%foldername = "app";
opufile = "opu" + filename;
cpufile = "cpu" + filename;
gpufile = "gpu" + filename;

stropu = "template <typename T> void " + opufile;
strgpu = "template <typename T>  __global__  void kernel" + gpufile;

tmp = "(T *f, T *xdg, T *udg1, T *udg2,  T *odg1, T *odg2,  T *wdg1, T *wdg2,  T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";

stropu = stropu + tmp + "{\n";
stropu = stropu + "\tfor (int i = 0; i <ng; i++) {\n";

strgpu = strgpu + tmp + "{\n";
strgpu = strgpu + "\tint i = threadIdx.x + blockIdx.x * blockDim.x;\n";
strgpu = strgpu + "\twhile (i<ng) {\n";

fstr = string(f(:));
str = "";
str = varsassign(str, "param", length(param), 0, fstr);
str = varsassign(str, "uinf", length(uinf), 0, fstr);
str = varsassign(str, "tau", length(tau), 0, fstr);
str = varsassign(str, "xdg", length(xdg), 1, fstr);
str = varsassign(str, "udg1", length(udg1), 1, fstr);
str = varsassign(str, "udg2", length(udg2), 1, fstr);
str = varsassign(str, "uhg", length(uhg), 1, fstr);
str = varsassign(str, "odg1", length(odg1), 1, fstr);
str = varsassign(str, "odg2", length(odg2), 1, fstr);
str = varsassign(str, "wdg1", length(wdg1), 1, fstr);
str = varsassign(str, "wdg2", length(wdg2), 1, fstr);
str = varsassign(str, "nlg", length(nlg), 1, fstr);
str = symsassign(str, f);

for i = 1:length(f(:))
    a = "f[" + num2str(i-1) + "*ng+i] =";
    b = "f[" + num2str(i-1) + "*ng+i] +=";
    str = strrep(str, a, b);
end

stropu = stropu + str + "\t}\n" + "}\n\n";
tmp = "template void " + opufile;
tmp = tmp + "(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int);\n";
tmp = tmp + "template void " + opufile;
tmp = tmp + "(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int);\n";
stropu = stropu + tmp;

fid = fopen(foldername + "/" + opufile + ".cpp", 'w');
fprintf(fid, char(stropu));
fclose(fid);

strgpu = strgpu + str + "\t\ti += blockDim.x * gridDim.x;\n";
strgpu = strgpu + "\t}\n" + "}\n\n";
tmp = "template <typename T> void " + gpufile;
tmp = tmp + "(T *f, T *xdg, T *udg1, T *udg2,  T *odg1, T *odg2,  T *wdg1, T *wdg2, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";
tmp = tmp + "{\n";
tmp = tmp + "\tint blockDim = 256;\n";
tmp = tmp + "\tint gridDim = (ng + blockDim - 1) / blockDim;\n";
tmp = tmp + "\tgridDim = (gridDim>1024)? 1024 : gridDim;\n";
tmp = tmp + "\tkernel" + gpufile + "<<<gridDim, blockDim>>>(f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);\n";
tmp = tmp + "}\n\n";
tmp = tmp + "template void " + gpufile;
tmp = tmp + "(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double time, int, int, int, int, int, int, int, int);\n";
tmp = tmp + "template void " + gpufile;
tmp = tmp + "(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float time, int, int, int, int, int, int, int, int);\n";
strgpu = strgpu + tmp;

fid = fopen(foldername + "/" + gpufile + ".cu", 'w');
fprintf(fid, char(strgpu));
fclose(fid);

strcpu = strrep(stropu, 'opu', "cpu");
strcpu = strrep(strcpu, "for (int i = 0; i <ng; i++) {", "#pragma omp parallel for\n\tfor (int i = 0; i <ng; i++) {");

iocpu = fopen(foldername + "/" + cpufile + ".cpp", 'w');
fprintf(fid, char(strcpu));
fclose(iocpu);

end
