import os

def nocodeelem2(filename, foldername):
    cpufile = "Kokkos" + filename
    tmp = "(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw, const int nce, const int npe, const int ne)\n"
    strkk = "void " + cpufile
    strkk += tmp + "{\n"
    strkk += "}\n"

    with open(os.path.join(foldername, cpufile + ".cpp"), "w") as fid:
        fid.write(strkk)
        
# def nocodeelem2(filename):

#     foldername = "app";
#     opufile = "opu" + filename;
#     cpufile = "cpu" + filename;
#     gpufile = "gpu" + filename;

#     stropu = "template <typename T> void " + opufile;
    
#     tmp = "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne)\n";

#     stropu = stropu + tmp + "{\n";
#     stropu = stropu + "}\n\n";

#     tmp = "template void " + opufile;
#     tmp = tmp + "(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int, int, int);\n";
#     tmp = tmp + "template void " + opufile;
#     tmp = tmp + "(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int, int, int);\n";

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
