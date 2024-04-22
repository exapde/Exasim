function nocodeelem2(filename::String)

foldername = "app";
opufile = "opu" * filename;
cpufile = "cpu" * filename;
gpufile = "gpu" * filename;

stropu = "template <typename T> void " * opufile;

tmp = "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne)\n";

stropu = stropu * tmp * "{\n";
stropu = stropu * "}\n\n";

tmp = "template void " * opufile;
tmp = tmp * "(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int, int, int);\n";
tmp = tmp * "template void " * opufile;
tmp = tmp * "(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int, int, int);\n";

stropu = stropu * tmp;
strcpu = replace(stropu, "opu" => "cpu");
strgpu = replace(stropu, "opu" => "gpu");

strgpu = strgpu * "\ntemplate <typename T> void __device__ device" * gpufile;

tmp = "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne)\n";

tmp =  tmp * "{\n";
tmp = tmp * "}\n\n";
strgpu = strgpu * tmp;

ioopu = open(foldername * "/" * opufile * ".cpp", "w");
write(ioopu, stropu);
close(ioopu);

iocpu = open(foldername * "/" * cpufile * ".cpp", "w");
write(iocpu, strcpu);
close(iocpu);

iogpu = open(foldername * "/" * gpufile * ".cu", "w");
write(iogpu, strgpu);
close(iogpu);

end
