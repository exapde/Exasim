from varsassign import varsassign
from sympyassign import sympyassign

def gencodebou(filename, f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time):

    opufile = "opu" + filename;
    gpufile = "gpu" + filename;

    stropu = "template <typename T> void " + opufile;
    strgpu = "template <typename T>  __global__  void kernel" + gpufile;
    
    tmp = "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";

    stropu = stropu + tmp + "{\n";
    stropu = stropu +  "\tfor (int i = 0; i <ng; i++) {\n";

    strgpu = strgpu + tmp + "{\n";
    strgpu = strgpu + "\tint i = threadIdx.x + blockIdx.x * blockDim.x;\n";
    strgpu = strgpu + "\twhile (i<ng) {\n";

    str = "";
    str = varsassign(str, "param", len(param), 0);
    str = varsassign(str, "uinf", len(uinf), 0);
    str = varsassign(str, "tau", len(tau), 0);
    str = varsassign(str, "xdg", len(xdg), 1);
    str = varsassign(str, "udg", len(udg), 1);
    str = varsassign(str, "uhg", len(uhg), 1);
    str = varsassign(str, "odg", len(odg), 1);
    str = varsassign(str, "wdg", len(wdg), 1);
    str = varsassign(str, "nlg", len(nlg), 1);
    str = sympyassign(str, f);
    stropu = stropu +  str + "\t}\n" + "}\n\n";

    strgpu = strgpu + str + "\t\ti += blockDim.x * gridDim.x;\n";
    strgpu = strgpu + "\t}\n" + "}\n\n";

    return stropu, strgpu
