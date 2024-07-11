function gencodeelem(filename::String, f, xdg, udg, odg, wdg, uinf, param, time, foldername)

    cpufile = "Kokkos" * filename
    tmp = "(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n"
    str = "\tKokkos::parallel_for(" * "\"" * filename  * "\"" * ", ng, KOKKOS_LAMBDA(const size_t i) {\n"
    strkk = "void " * cpufile
    strkk = strkk * tmp * "{\n"
    strkk = strkk * str

    fstr = string(f[:])
    str = ""
    str = varsassign(str, "param", length(param), 0, fstr)
    str = varsassign(str, "uinf", length(uinf), 0, fstr)
    str = varsassign(str, "xdg", length(xdg), 1, fstr)
    str = varsassign(str, "udg", length(udg), 1, fstr)
    str = varsassign(str, "odg", length(odg), 1, fstr)
    str = varsassign(str, "wdg", length(wdg), 1, fstr)
    str = sympyassign(str, f)

    strkk = strkk * str * "\t});\n" * "}\n\n"
    strkk = replace(strkk, "T " => "dstype ")

    open(foldername * "/" * cpufile * ".cpp", "w") do fid
        write(fid, strkk)
    end

  end
