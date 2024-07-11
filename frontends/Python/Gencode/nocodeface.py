import os

def nocodeface(filename, foldername):
    cpufile = "Kokkos" + filename
    tmp = "(const dstype* f, const dstype* xdg, const dstype* udg1, const dstype* udg2, const dstype* odg1, const dstype* odg2, const dstype* wdg1, const dstype* wdg2, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n"
    strkk = "void " + cpufile
    strkk += tmp + "{\n"
    strkk += "}\n"

    with open(os.path.join(foldername, cpufile + ".cpp"), "w") as fid:
        fid.write(strkk)

# def nocodeface(filename):

#     foldername = "app";
#     opufile = "opu" + filename;
#     cpufile = "cpu" + filename;
#     gpufile = "gpu" + filename;

#     stropu = "template <typename T> void " + opufile;
    
#     tmp = "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";

#     stropu = stropu + tmp + "{\n";
#     stropu = stropu + "}\n\n";

#     tmp = "template void " + opufile;
#     tmp = tmp + "(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int);\n";
#     tmp = tmp + "template void " + opufile;
#     tmp = tmp + "(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int);\n";

#     stropu = stropu + tmp;
#     strcpu = stropu.replace("opu", "cpu");
#     strgpu = stropu.replace("opu", "gpu");

#     ioopu = open(foldername + "/" + opufile + ".cpp", "w");
#     ioopu.write(stropu);
#     ioopu.close();

#     iogpu = open(foldername + "/" + gpufile + ".cu", "w");
#     iogpu.write(strgpu);
#     iogpu.close();

#     iocpu = open(foldername + "/" + cpufile + ".cpp", "w");
#     iocpu.write(strcpu);
#     iocpu.close();

#     return 0
