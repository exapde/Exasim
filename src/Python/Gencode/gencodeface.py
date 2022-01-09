from gencodebou import gencodebou

def gencodeface(filename, f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time):

    foldername = "app";
    opufile = "opu" + filename;
    cpufile = "cpu" + filename;
    gpufile = "gpu" + filename;

    nbc = f.shape[1];

    stropu = "";
    strgpu = "";
    for k in range(1,nbc+1):
        str1, str2 = gencodebou(filename + str(k), f[:,k-1], xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time);
        stropu = stropu + str1;
        strgpu = strgpu + str2;

    tmp = "template <typename T> void " + opufile;
    tmp = tmp + "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";
    tmp = tmp + "{\n";
    for k in range(1,nbc+1):
        if k == 1:
            tmp = tmp + "\tif (ib == " + str(k) + ")\n";
        else:
            tmp = tmp + "\telse if (ib == " + str(k) + ")\n";
        tmp = tmp + "\t\t" + opufile + str(k) + "(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);\n";
    tmp = tmp + "}\n\n";

    tmp = tmp + "template void " + opufile;
    tmp = tmp + "(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int);\n";
    tmp = tmp + "template void " + opufile;
    tmp = tmp + "(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int);\n";
    stropu = stropu + tmp;

    tmp = "template <typename T> void " + gpufile;
    tmp = tmp + "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";
    tmp = tmp + "{\n";
    tmp = tmp + "\tint blockDim = 256;\n";
    tmp = tmp + "\tint gridDim = (ng + blockDim - 1) / blockDim;\n";
    tmp = tmp + "\tgridDim = (gridDim>1024)? 1024 : gridDim;\n";
    for k in range(1,nbc+1):
        if k == 1:
            tmp = tmp + "\tif (ib == " + str(k) + ")\n";
        elif k<nbc:
            tmp = tmp + "\telse if (ib == " + str(k) + ")\n";
        else:
            tmp = tmp + "\telse if (ib == " + str(k) + ")\n";
        tmp = tmp + "\t\tkernel"  + gpufile + str(k) + "<<<gridDim, blockDim>>>(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);\n";
    tmp = tmp + "}\n\n";

    tmp = tmp + "template void " + gpufile;
    tmp = tmp + "(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int);\n";
    tmp = tmp + "template void " + gpufile;
    tmp = tmp + "(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int);\n";
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
