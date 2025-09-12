function [stropu, strgpu] = gencodebou(filename, f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time)

opufile = "opu" + filename;
gpufile = "gpu" + filename;

stropu = "template <typename T> void " + opufile;
strgpu = "template <typename T>  __device__  void device" + gpufile;

tmp = "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";

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
str = varsassign(str, "udg", length(udg), 1, fstr);
str = varsassign(str, "uhg", length(uhg), 1, fstr);
str = varsassign(str, "odg", length(odg), 1, fstr);
str = varsassign(str, "wdg", length(wdg), 1, fstr);
str = varsassign(str, "nlg", length(nlg), 1, fstr);
str = symsassign(str, f);
stropu = stropu + str + "\t}\n" + "}\n\n";

strgpu = strgpu + str + "\t\ti += blockDim.x * gridDim.x;\n";
strgpu = strgpu + "\t}\n" + "}\n\n";

strgpu = strgpu + "template <typename T>  __global__  void kernel" + gpufile;
tmp = "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";
tmp = tmp + "{\n";
tmp = tmp + "\tdevice" + gpufile + "(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);\n";
tmp = tmp + "}\n\n";
strgpu = strgpu + tmp;

end



