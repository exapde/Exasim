import sympy, os
from varsassign import varsassign

def gencodeelem4(filename, f, xdg, uinf, param, foldername):
    
    cpufile = "cpu" + filename
    tmp = "(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)\n"
    str_code = "\tfor (int i = 0; i < ng; i++) {\n"
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

    n = len(f)
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

    strkk += str_code + "\t}\n" + "}\n\n"
    strkk = strkk.replace("T ", "dstype ")

    with open(os.path.join(foldername, cpufile + ".cpp"), "w") as fid:
        fid.write(strkk)

