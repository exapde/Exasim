import sympy
from varsassign import varsassign

def gencodeelem2(filename, f, xdg, udg, odg, wdg, uinf, param, time, foldername):
    cpufile = "Kokkos" + filename
    tmp = "(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw, const int nce, const int npe, const int ne)\n"
    str_code = "\tKokkos::parallel_for(" + "\"" + filename + "\"" + ", ng, KOKKOS_LAMBDA(const size_t i) {\n"
    strkk = "void " + cpufile
    strkk += tmp + "{\n"
    strkk += str_code

    str_code = "\t\tint j = i % npe;\n"
    str_code += "\t\tint k = i / npe;\n"

    fstr = str(f.flatten())
    str_code = varsassign(str_code, "param", len(param), 0, fstr)
    str_code = varsassign(str_code, "uinf", len(uinf), 0, fstr)

    varname = "xdg"
    for i in range(1, len(xdg) + 1):
        str1 = varname + str(i)
        if str1 in fstr:
            str2 = varname + "[j+npe*" + str(i - 1) + "+npe*ncx*k]"
            str_code += "\t\tdstype {} = {};\n".format(str1, str2)

    varname = "udg"
    for i in range(1, len(udg) + 1):
        str1 = varname + str(i)
        if str1 in fstr:
            str2 = varname + "[j+npe*" + str(i - 1) + "+npe*nc*k]"
            str_code += "\t\tdstype {} = {};\n".format(str1, str2)

    varname = "odg"
    for i in range(1, len(odg) + 1):
        str1 = varname + str(i)
        if str1 in fstr:
            str2 = varname + "[j+npe*" + str(i - 1) + "+npe*nco*k]"
            str_code += "\t\tdstype {} = {};\n".format(str1, str2)

    varname = "wdg"
    for i in range(1, len(wdg) + 1):
        str1 = varname + str(i)
        if str1 in fstr:
            str2 = varname + "[j+npe*" + str(i - 1) + "+npe*ncw*k]"
            str_code += "\t\tdstype {} = {};\n".format(str1, str2)

    n = len(f)
    # Simplify f if needed using sympy
    ts, fs = sympy.cse(f)
    m = len(ts)
    for i in range(m):
        str1 = sympy.ccode(ts[i][0])
        str2 = sympy.ccode(ts[i][1])
        str_code += "\t\tdstype {} = {};\n".format(str1, str2)

    for i in range(n):
        str1 = "\t\tf[j+npe*" + str(i) + "+npe*nce*k]"
        str2 = sympy.ccode(fs[i])
        str_code += "{} = {};\n".format(str1, str2)

    strkk += str_code + "\t});\n" + "}\n\n"
    strkk = strkk.replace("T ", "dstype ")

    with open(os.path.join(foldername, cpufile + ".cpp"), "w") as fid:
        fid.write(strkk)

    return 0
    
# def gencodeelem2(filename, f, xdg, udg, odg, wdg, uinf, param, time):

#     foldername = "app";
#     opufile = "opu" + filename;
#     cpufile = "cpu" + filename;
#     gpufile = "gpu" + filename;

#     stropu = "template <typename T> void " + opufile;
#     strgpu = "template <typename T>  __global__  void kernel" + gpufile;
    
#     tmp = "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne)\n";

#     stropu = stropu + tmp + "{\n";
#     stropu = stropu + "\tfor (int i = 0; i <ng; i++) {\n";

#     strgpu = strgpu + tmp + "{\n";
#     strgpu = strgpu + "\tint i = threadIdx.x + blockIdx.x * blockDim.x;\n";
#     strgpu = strgpu + "\twhile (i<ng) {\n";

#     mystr = "\t\tint j = i%npe;\n";
#     mystr = mystr + "\t\tint k = (i-j)/npe;\n";

#     mystr = varsassign(mystr, "param", len(param), 0);
#     mystr = varsassign(mystr, "uinf", len(uinf), 0);

#     varname = "xdg";
#     for i in range(0, len(xdg)):
#         str1 = varname + str(i+1);
#         str2 = varname + "[j+npe*" + str(i) + "+npe*ncx*k" + "]";
#         mystr = mystr + "\t\tT " + str1 + " = " + str2 + ";\n";

#     varname = "udg";
#     for i in range(0, len(udg)):
#         str1 = varname + str(i+1);
#         str2 = varname + "[j+npe*" + str(i) + "+npe*nc*k" + "]";
#         mystr = mystr + "\t\tT " + str1 + " = " + str2 + ";\n";

#     varname = "odg";
#     for i in range(0, len(odg)):
#         str1 = varname + str(i+1);
#         str2 = varname + "[j+npe*" + str(i) + "+npe*nco*k" + "]";
#         mystr = mystr + "\t\tT " + str1 + " = " + str2 + ";\n";

#     varname = "wdg";
#     for i in range(0, len(wdg)):
#         str1 = varname + str(i+1);
#         str2 = varname + "[j+npe*" + str(i) + "+npe*ncw*k" + "]";
#         mystr = mystr + "\t\tT " + str1 + " = " + str2 + ";\n";

#     n = len(f);
#     #for i in range(0,n):
#     #    f[i] = sympy.simplify(f[i]);

#     ts, fs = sympy.cse(f);
#     m = len(ts);
#     for i in range(0, m):
#         str1 = sympy.ccode(ts[i][0]);
#         str2 = sympy.ccode(ts[i][1]);
#         mystr = mystr + "\t\tT " + str1 + " = " + str2 + ";\n";

#     for i in range(0, n):
#         str1 = "\t\tf[j+npe*" + str(i) + "+npe*nce*k" + "]";
#         str2 = sympy.ccode(fs[i]);
#         mystr = mystr + str1 + " = " + str2 + ";\n";

#     stropu = stropu + mystr + "\t}\n" + "}\n\n";
#     tmp = "template void " + opufile;
#     tmp = tmp + "(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int, int, int);\n";
#     tmp = tmp + "template void " + opufile;
#     tmp = tmp + "(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int, int, int);\n";
#     stropu = stropu + tmp;

#     strgpu = strgpu + mystr + "\t\ti += blockDim.x * gridDim.x;\n";
#     strgpu = strgpu + "\t}\n" + "}\n\n";
#     tmp = "template <typename T> void " + gpufile;
#     tmp = tmp + "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne)\n";
#     tmp = tmp + "{\n";
#     tmp = tmp + "\tint blockDim = 256;\n";
#     tmp = tmp + "\tint gridDim = (ng + blockDim - 1) / blockDim;\n";
#     tmp = tmp + "\tgridDim = (gridDim>1024)? 1024 : gridDim;\n";
#     tmp = tmp + "\tkernel" + gpufile + "<<<gridDim, blockDim>>>(f, xdg, udg, odg, wdg, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne);\n";
#     tmp = tmp + "}\n\n";
#     tmp = tmp + "template void " + gpufile;
#     tmp = tmp + "(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int, int, int);\n";
#     tmp = tmp + "template void " + gpufile;
#     tmp = tmp + "(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int, int, int);\n";
#     strgpu = strgpu + tmp;

#     ioopu = open(foldername + "/" + opufile + ".cpp", "w");
#     ioopu.write(stropu);
#     ioopu.close();

#     iogpu = open(foldername + "/" + gpufile + ".cu", "w");
#     iogpu.write(strgpu);
#     iogpu.close();

#     strcpu = stropu.replace("opu", "cpu");
#     strcpu.replace("for (int i = 0; i <ng; i++) {", "#pragma omp parallel for\n\tfor (int i = 0; i <ng; i++) {");

#     iocpu = open(foldername + "/" + cpufile + ".cpp", "w");
#     iocpu.write(strcpu);
#     iocpu.close();

#     return 0
