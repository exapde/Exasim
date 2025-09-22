function hdgnocodeface2(filename::String, foldername::String)

    cpufile = "Hdg" * filename
    tmp = "(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n"
    strkk = "void " * cpufile
    strkk = strkk * tmp * "{\n"
    strkk = strkk * "}\n"

    open(foldername * "/" * cpufile * ".cpp", "w") do fid
        write(fid, strkk)
    end
end