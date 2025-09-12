function gencodeface3(filename::String, f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, foldername)

cpufile = "Kokkos" * filename
tmp = "(dstype* f, const dstype* xdg, const dstype* udg1, const dstype* udg2, const dstype* odg1, const dstype* odg2, const dstype* wdg1, const dstype* wdg2, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n"
str = "\tKokkos::parallel_for(" * "\"" * filename * "\"" * ", ng, KOKKOS_LAMBDA(const size_t i) {\n"
strkk = "void " * cpufile
strkk = strkk * tmp * "{\n"
strkk = strkk * str

fstr = string(f[:])
str = ""
str = varsassign(str, "param", length(param), 0, fstr)
str = varsassign(str, "uinf", length(uinf), 0, fstr)
str = varsassign(str, "tau", length(tau), 0, fstr)
str = varsassign(str, "xdg", length(xdg), 1, fstr)
str = varsassign(str, "udg1", length(udg1), 1, fstr)
str = varsassign(str, "udg2", length(udg2), 1, fstr)
str = varsassign(str, "uhg", length(uhg), 1, fstr)
str = varsassign(str, "odg1", length(odg1), 1, fstr)
str = varsassign(str, "odg2", length(odg2), 1, fstr)
str = varsassign(str, "wdg1", length(wdg1), 1, fstr)
str = varsassign(str, "wdg2", length(wdg2), 1, fstr)
str = varsassign(str, "nlg", length(nlg), 1, fstr)
str = sympyassign(str, f)

for i in 1:length(f[:])
  a = "f[" * string(i-1) * "*ng+i] ="
  b = "f[" * string(i-1) * "*ng+i] +="
  str = replace(str, a => b)
end

strkk = strkk * str * "\t});\n" * "}\n\n"
strkk = replace(strkk, "T ", "dstype ")

open(foldername * "/" * cpufile * ".cpp", "w") do fid
    write(fid, strkk)
end

end
