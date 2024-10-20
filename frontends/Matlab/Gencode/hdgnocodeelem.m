function strkk = hdgnocodeelem(filename, foldername)

cpufile = "Hdg" + filename;
tmp = "(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n";
strkk = "void " + cpufile;
strkk = strkk + tmp + "{\n";

strkk = strkk + "}\n"; 
fid = fopen(foldername + "/" + cpufile + ".cpp", 'w');
fprintf(fid, char(strkk));
fclose(fid);

end

