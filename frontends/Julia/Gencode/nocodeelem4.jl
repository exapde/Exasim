function nocodeelem4(filename::String, foldername::String)

cpufile = "cpu" * filename
tmp = "(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)\n"
strkk = "void " * cpufile
strkk = strkk * tmp * "{\n"
strkk = strkk * "}\n"

open(foldername * "/" * cpufile * ".cpp", "w") do fid
    write(fid, strkk)
end

end
