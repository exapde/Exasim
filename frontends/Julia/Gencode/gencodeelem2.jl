function gencodeelem2(filename::String, f, xdg, udg, odg, wdg, uinf, param, time, foldername)

  cpufile = "Kokkos" * filename
  tmp = "(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw, const int nce, const int npe, const int ne)\n"
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

  varname = "udg"
  for i in 1:length(udg)
      str1 = varname * string(i)
      if occursin(str1, fstr)
          str2 = varname * "[j+npe*" * string(i-1) * "+npe*nc*k]"
          str = str * "\t\tdstype " * str1 * " = " * str2 * ";\n"
      end
  end

  varname = "odg"
  for i in 1:length(odg)
      str1 = varname * string(i)
      if occursin(str1, fstr)
          str2 = varname * "[j+npe*" * string(i-1) * "+npe*nco*k]"
          str = str * "\t\tdstype " * str1 * " = " * str2 * ";\n"
      end
  end

  varname = "wdg"
  for i in 1:length(wdg)
      str1 = varname * string(i)
      if occursin(str1, fstr)
          str2 = varname * "[j+npe*" * string(i-1) * "+npe*ncw*k]"
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
  for i = 1:n
      str1 = "\t\tf[j+npe*" * string(i-1) * "+npe*nce*k" * "]";
      str2 = string(sympy.ccode(fs[1][i]));
      str = str * str1 * " = " * str2 * ";\n";
  end
  
  strkk = strkk * str * "\t});\n" * "}\n\n"
  strkk = replace(strkk, "dstype " => "dstype ")
  open(foldername * "/" * cpufile * ".cpp", "w") do fid
      write(fid, strkk)
  end

end

