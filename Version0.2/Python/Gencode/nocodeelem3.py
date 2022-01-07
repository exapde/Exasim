def nocodeelem3(filename):

    foldername = "app";
    opufile = "opu" + filename;
    cpufile = "cpu" + filename;
    gpufile = "gpu" + filename;

    stropu = "template <typename T> void " + opufile;
    
    tmp = "(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)\n";

    stropu = stropu + tmp + "{\n";
    stropu = stropu + "}\n\n";

    tmp = "template void " + opufile;
    tmp = tmp + "(double *, double *, double *, double *, int, int, int, int, int, int);\n";
    tmp = tmp + "template void " + opufile;
    tmp = tmp + "(float *, float *, float *, float *, int, int, int, int, int, int);\n";

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
