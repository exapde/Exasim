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

    # Enzyme
    tmp = "#ifdef _ENZYME\n"
    for k in range(1,nbc+1):
        if gpufile == "gpuUbou":
            tmp = tmp + "template <typename T> __global__ void kernelGrad" + gpufile  + str(k) + "Enzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *wdg, T *dwdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";
        else:
            tmp = tmp + "template <typename T> __global__ void kernelGrad" + gpufile + str(k) + "Enzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *dodg, T *wdg, T *dwdg, T *uhg, T *duhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";
        
        tmp = tmp + "{\n";
        tmp = tmp + "\t__enzyme_fwddiff" + gpufile + "((void*)device"  + gpufile + str(k) + "<T>,\n"
        tmp = tmp + "\t\t\t  enzyme_dup, f, df,\n"
        tmp = tmp + "\t\t\t enzyme_const, xg,\n"
        tmp = tmp + "\t\t\t enzyme_dup, udg, dudg,\n"
        if gpufile == "gpuUbou":
            tmp = tmp + "\t\t\t enzyme_const, odg,\n"
        else:
            tmp = tmp + "\t\t\t enzyme_dup, odg, dodg,\n"
    
        tmp = tmp + "\t\t\t enzyme_dup, wdg, dwdg,\n"
        if gpufile == "gpuUbou":
            tmp = tmp + "\t\t\t enzyme_const, uhg,\n"
        else:
            tmp = tmp + "\t\t\t enzyme_dup, uhg, duhg,\n"

        tmp = tmp + "\t\t\t enzyme_const, nlg,\n"
        tmp = tmp + "\t\t\t enzyme_const, tau,\n"
        tmp = tmp + "\t\t\t enzyme_const, uinf,\n"
        tmp = tmp + "\t\t\t enzyme_const, param,\n"
        tmp = tmp + "\t\t\t enzyme_const, time,\n"
        tmp = tmp + "\t\t\t enzyme_const, modelnumber,\n"
        tmp = tmp + "\t\t\t enzyme_const, ng,\n"
        tmp = tmp + "\t\t\t enzyme_const, nc,\n"
        tmp = tmp + "\t\t\t enzyme_const, ncu,\n"
        tmp = tmp + "\t\t\t enzyme_const, nd,\n"
        tmp = tmp + "\t\t\t enzyme_const, ncx,\n"
        tmp = tmp + "\t\t\t enzyme_const, nco,\n"
        tmp = tmp + "\t\t\t enzyme_const, ncw);\n"
    tmp = tmp + "}\n\n";

    tmp = tmp + "template <typename T> void " + gpufile + "Enzyme";
    if gpufile == "gpuUbou":
        tmp = tmp + "(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *wdg, T *dwdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";
    else:
        tmp = tmp + "(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *dodg, T *wdg, T *dwdg, T *uhg, T *duhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";

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

        if gpufile == "gpuUbou":
            tmp = tmp + "\t\tkernelGrad"  + gpufile + str(k) + "Enzyme<<<gridDim, blockDim>>>(f, df, xg, udg, dudg, odg, wdg, dwdg, uhg, nlg, tau, uinf, param, time, modelnumber, ib, ng, nc, ncu, nd, ncx, nco, ncw);\n";
        else:
            tmp = tmp + "\t\tkernelGrad"  + gpufile + str(k) + "Enzyme<<<gridDim, blockDim>>>(f, df, xg, udg, dudg, odg, dodg, wdg, dwdg, uhg, duhg, nlg, tau, uinf, param, time, modelnumber, ib, ng, nc, ncu, nd, ncx, nco, ncw);\n";

    tmp = tmp + "}\n\n";

    tmp = tmp + "template void " + gpufile + "Enzyme";
    if gpufile == "gpuUbou":
        tmp = tmp + "(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int);\n";
    else:
        tmp = tmp + "(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int);\n";

    tmp = tmp + "#endif"
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
