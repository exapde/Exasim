from varsassign import varsassign
from sympyassign import sympyassign

def gencodebou(filename, f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time):
    cpufile = "Kokkos" + filename
    tmp = "(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n"
    str_code = "\tKokkos::parallel_for(" + "\"" + filename + "\"" + ", ng, KOKKOS_LAMBDA(const size_t i) {\n"
    strkk = "void " + cpufile
    strkk += tmp + "{\n"
    strkk += str_code

    fstr = str(f.flatten())
    str_code = ""
    str_code = varsassign(str_code, "param", len(param), 0, fstr)
    str_code = varsassign(str_code, "uinf", len(uinf), 0, fstr)
    str_code = varsassign(str_code, "tau", len(tau), 0, fstr)
    str_code = varsassign(str_code, "xdg", len(xdg), 1, fstr)
    str_code = varsassign(str_code, "udg", len(udg), 1, fstr)
    str_code = varsassign(str_code, "uhg", len(uhg), 1, fstr)
    str_code = varsassign(str_code, "odg", len(odg), 1, fstr)
    str_code = varsassign(str_code, "wdg", len(wdg), 1, fstr)
    str_code = varsassign(str_code, "nlg", len(nlg), 1, fstr)
    str_code = sympyassign(str_code, f)

    strkk += str_code + "\t});\n" + "}\n\n"
    strkk = strkk.replace("T ", "dstype ")

    return strkk

# def gencodebou(filename, f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time):

#     opufile = "opu" + filename;
#     gpufile = "gpu" + filename;

#     stropu = "template <typename T> void " + opufile;
#     strgpu = "template <typename T>  __device__  void device" + gpufile;
    
#     tmp = "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";

#     stropu = stropu + tmp + "{\n";
#     stropu = stropu +  "\tfor (int i = 0; i <ng; i++) {\n";

#     strgpu = strgpu + tmp + "{\n";
#     strgpu = strgpu + "\tint i = threadIdx.x + blockIdx.x * blockDim.x;\n";
#     strgpu = strgpu + "\twhile (i<ng) {\n";

#     str = "";
#     str = varsassign(str, "param", len(param), 0);
#     str = varsassign(str, "uinf", len(uinf), 0);
#     str = varsassign(str, "tau", len(tau), 0);
#     str = varsassign(str, "xdg", len(xdg), 1);
#     str = varsassign(str, "udg", len(udg), 1);
#     str = varsassign(str, "uhg", len(uhg), 1);
#     str = varsassign(str, "odg", len(odg), 1);
#     str = varsassign(str, "wdg", len(wdg), 1);
#     str = varsassign(str, "nlg", len(nlg), 1);
#     str = sympyassign(str, f);
#     stropu = stropu +  str + "\t}\n" + "}\n\n";

#     strgpu = strgpu + str + "\t\ti += blockDim.x * gridDim.x;\n";
#     strgpu = strgpu + "\t}\n" + "}\n\n";

#     strgpu = strgpu + "template <typename T>  __global__  void kernel" + gpufile
#     tmp = "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";
#     tmp = tmp + "{\n";
#     tmp = tmp + "\tdevice" + gpufile + "(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);\n";
#     tmp = tmp + "}\n\n";
#     strgpu = strgpu + tmp
    
#     return stropu, strgpu
