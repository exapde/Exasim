function nocodeelem3(filename)

foldername = "app";
opufile = "opu" + filename;
cpufile = "cpu" + filename;
gpufile = "gpu" + filename;

stropu = "template <typename T> void " + opufile;

tmp = "(T *f, T *xdg, T *uinf, T *param, int ng, int ncx, int nce, int npe, int ne)\n";

stropu = stropu + tmp + "{\n";
stropu = stropu + "}\n\n";

tmp = "template void " + opufile;
tmp = tmp + "(double *, double *, double *, double *, int, int, int, int, int);\n";
tmp = tmp + "template void " + opufile;
tmp = tmp + "(float *, float *, float *, float *, int, int, int, int, int);\n";

stropu = stropu + tmp;
strcpu = strrep(stropu, "opu", "cpu");
strgpu = strrep(stropu, "opu", "gpu");

fid = fopen(foldername + "/" + opufile + ".cpp", 'w');
fprintf(fid, char(stropu));
fclose(fid);

fid = fopen(foldername + "/" + gpufile + ".cu", 'w');
fprintf(fid, char(strgpu));
fclose(fid);

iocpu = fopen(foldername + "/" + cpufile + ".cpp", 'w');
fprintf(fid, char(strcpu));
fclose(iocpu);

end
