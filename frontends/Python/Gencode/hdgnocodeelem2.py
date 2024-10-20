def hdgnocodeelem2(filename, foldername):
    cpufile = "Hdg" + filename
    tmp = "(dstype* f, dstype* f_wdg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n"
    strkk = "void " + cpufile
    strkk += tmp + "{\n"
    strkk += "}\n"

    with open(foldername + "/" + cpufile + ".cpp", "w") as fid:
        fid.write(strkk)

# Example usage:
# foldername = "example_folder"
# filename = "example_filename"
# hdgnocodeelem2(filename, foldername)
