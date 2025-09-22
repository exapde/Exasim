from hdggencodebou import hdggencodebou

def hdggencodeface(filename, f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, foldername):
    strkk = ""
    nbc = f.shape[1]
    for k in range(nbc):
        str1 = hdggencodebou(filename + str(k + 1), f[:, k], xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time)
        strkk += str1

    cpufile = "Hdg" + filename
    tmp = "(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n"
    tmp = "void " + cpufile + tmp
    tmp += "{\n"
    for k in range(nbc):
        if k == 0:
            tmp += "\tif (ib == " + str(k + 1) + ")\n"
        else:
            tmp += "\telse if (ib == " + str(k + 1) + ")\n"
        tmp += "\t\t" + cpufile + str(k + 1) + "(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);\n"
    tmp += "}\n\n"

    strkk += tmp
    with open(foldername + "/" + cpufile + ".cpp", "w") as fid:
        fid.write(strkk)

