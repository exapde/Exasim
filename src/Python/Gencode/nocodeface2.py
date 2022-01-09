def nocodeface2(filename):

    foldername = "app";
    opufile = "opu" + filename;
    cpufile = "cpu" + filename;
    gpufile = "gpu" + filename;

    stropu = "template <typename T> void " + opufile;
    
    tmp = "(T *f, T *xdg, T *udg1, T *udg2, T *odg1, T *odg2,  T *wdg1, T *wdg2, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";

    stropu = stropu + tmp + "{\n";
    stropu = stropu + "}\n\n";

    tmp = "template void " + opufile;
    tmp = tmp + "(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int);\n";
    tmp = tmp + "template void " + opufile;
    tmp = tmp + "(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int);\n";

    stropu = stropu + tmp;
    strcpu = stropu.replace("opu", "cpu");
    strgpu = stropu.replace("opu", "gpu");

    ioopu = open(foldername + "/" + opufile + ".cpp", "w");
    ioopu.write(stropu);
    ioopu.close();

    iogpu = open(foldername + "/" + gpufile + ".cu", "w");
    iogpu.write(strgpu);
    iogpu.close();

    iocpu = open(foldername + "/" + cpufile + ".cpp", "w");
    iocpu.write(strcpu);
    iocpu.close();

    return 0
