function strkk = kknocodeface(filename, foldername)

cpufile = "Kokkos" + filename;
tmp = "(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n";
strkk = "void " + cpufile;
strkk = strkk + tmp + "{\n";

strkk = strkk + "}\n"; 
fid = fopen(foldername + "/" + cpufile + ".cpp", 'w');
fprintf(fid, char(strkk));
fclose(fid);

end
