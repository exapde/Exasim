import os

def nocodeelem4(filename, foldername):
    cpufile = "cpu" + filename
    tmp = "(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)\n"
    strkk = "void " + cpufile
    strkk += tmp + "{\n"
    strkk += "}\n"

    with open(os.path.join(foldername, cpufile + ".cpp"), "w") as fid:
        fid.write(strkk)