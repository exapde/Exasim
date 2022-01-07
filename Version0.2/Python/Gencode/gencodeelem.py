from varsassign import varsassign
from sympyassign import sympyassign

def gencodeelem(filename, f, xdg, udg, odg, wdg, uinf, param, time):

    foldername = "app";
    opufile = "opu" + filename;
    cpufile = "cpu" + filename;
    gpufile = "gpu" + filename;

    stropu = "template <typename T> void " + opufile;
    strgpu = "template <typename T>  __global__  void kernel" + gpufile;
    
    tmp = "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";

    stropu = stropu + tmp + "{\n";
    stropu = stropu + "\tfor (int i = 0; i <ng; i++) {\n";

    strgpu = strgpu + tmp + "{\n";
    strgpu = strgpu + "\tint i = threadIdx.x + blockIdx.x * blockDim.x;\n";
    strgpu = strgpu + "\twhile (i<ng) {\n";

    mystr = "";
    mystr = varsassign(mystr, "param", len(param), 0);
    mystr = varsassign(mystr, "uinf", len(uinf), 0);
    mystr = varsassign(mystr, "xdg", len(xdg), 1);
    mystr = varsassign(mystr, "udg", len(udg), 1);
    mystr = varsassign(mystr, "odg", len(odg), 1);
    mystr = varsassign(mystr, "wdg", len(wdg), 1);
    mystr = sympyassign(mystr, f);

    stropu = stropu + mystr + "\t}\n" + "}\n\n";
    tmp = "template void " + opufile;
    tmp = tmp + "(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int);\n";
    tmp = tmp + "template void " + opufile;
    tmp = tmp + "(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int);\n";
    stropu = stropu + tmp;

    strgpu = strgpu + mystr + "\t\ti += blockDim.x * gridDim.x;\n";
    strgpu = strgpu + "\t}\n" + "}\n\n";
    tmp = "template <typename T> void " + gpufile;
    tmp = tmp + "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";    
    tmp = tmp + "{\n";
    tmp = tmp + "\tint blockDim = 256;\n";
    tmp = tmp + "\tint gridDim = (ng + blockDim - 1) / blockDim;\n";
    tmp = tmp + "\tgridDim = (gridDim>1024)? 1024 : gridDim;\n";
    tmp = tmp + "\tkernel" + gpufile + "<<<gridDim, blockDim>>>(f, xdg, udg, odg, wdg, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);\n";
    tmp = tmp + "}\n\n";
    tmp = tmp + "template void " + gpufile;
    tmp = tmp + "(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int);\n";
    tmp = tmp + "template void " + gpufile;
    tmp = tmp + "(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int);\n";
    strgpu = strgpu + tmp;

    ioopu = open(foldername + "/" + opufile + ".cpp", "w");
    ioopu.write(stropu);
    ioopu.close();

    iogpu = open(foldername + "/" + gpufile + ".cu", "w");
    iogpu.write(strgpu);
    iogpu.close();

    strcpu = stropu.replace("opu", "cpu");
    strcpu.replace("for (int i = 0; i <ng; i++) {", "#pragma omp parallel for\n\tfor (int i = 0; i <ng; i++) {");

    iocpu = open(foldername + "/" + cpufile + ".cpp", "w");
    iocpu.write(strcpu);
    iocpu.close();

    return 0
