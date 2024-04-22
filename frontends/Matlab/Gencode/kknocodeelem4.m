function strkk = kknocodeelem4(filename, foldername)

cpufile = "cpu" + filename;
tmp = "(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)\n";
strkk = "void " + cpufile;
strkk = strkk + tmp + "{\n";

strkk = strkk + "}\n"; 
fid = fopen(foldername + "/" + cpufile + ".cpp", 'w');
fprintf(fid, char(strkk));
fclose(fid);

end
