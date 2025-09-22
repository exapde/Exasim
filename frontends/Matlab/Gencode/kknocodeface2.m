function strkk = kknocodeface2(filename, foldername)

cpufile = "Kokkos" + filename;
tmp = "(const dstype* f, const dstype* xdg, const dstype* udg1, const dstype* udg2,  const dstype* odg1, const dstype* odg2,  const dstype* wdg1, const dstype* wdg2,  const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n";
strkk = "void " + cpufile;
strkk = strkk + tmp + "{\n";

strkk = strkk + "}\n"; 
fid = fopen(foldername + "/" + cpufile + ".cpp", 'w');
fprintf(fid, char(strkk));
fclose(fid);

end
