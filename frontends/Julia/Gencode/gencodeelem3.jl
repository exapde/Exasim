function gencodeelem3(filename::String, f, xdg, uinf, param, foldername)

cpufile = "Kokkos" * filename
tmp = "(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)\n"
str = "\tKokkos::parallel_for(" * "\"" * filename * "\"" * ", ng, KOKKOS_LAMBDA(const size_t i) {\n"
strkk = "void " * cpufile
strkk = strkk * tmp * "{\n"
strkk = strkk * str

str = "\t\tint j = i % npe;\n"
str = str * "\t\tint k = i / npe;\n"

fstr = string(f[:])
str = varsassign(str, "param", length(param), 0, fstr)
str = varsassign(str, "uinf", length(uinf), 0, fstr)
varname = "xdg"
for i in 1:length(xdg)
    str1 = varname * string(i)
    if occursin(str1, fstr)
        str2 = varname * "[j+npe*" * string(i-1) * "+npe*ncx*k]"
        str = str * "\t\tdstype " * str1 * " = " * str2 * ";\n"
    end
end

n = length(f);
#for i=1:n
#    f[i] = sympy.simplify(f[i]);
#end
ts, fs = sympy.cse(f);
m = length(ts);
for i = 1:m
    str1 = string(sympy.ccode(ts[i][1]));
    str2 = string(sympy.ccode(ts[i][2]));
    str = str * "\t\tT " * str1 * " = " * str2 * ";\n";
end
if length(fs)==1
    for i = 1:n
        str1 = "\t\tf[j+npe*" * string(i-1) * "+npe*nce*k" * "]";
        str2 = string(sympy.ccode(fs[1][i]));
        str = str * str1 * " = " * str2 * ";\n";
    end
elseif length(fs)>1
    for i = 1:n
        str1 = "\t\tf[j+npe*" * string(i-1) * "+npe*nce*k" * "]";
        str2 = string(sympy.ccode(fs[i]));
        str = str * str1 * " = " * str2 * ";\n";
    end
else
    error("Symbolic expressions are wrong");
end

strkk = strkk * str * "\t});\n" * "}\n\n"
strkk = replace(strkk, "T " => "dstype ")
open(foldername * "/" * cpufile * ".cpp", "w") do fid
    write(fid, strkk)
end

end
