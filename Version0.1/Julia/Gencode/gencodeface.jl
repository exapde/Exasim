function gencodeface(filename::String, f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time)

foldername = "app";
opufile = "opu" * filename;
cpufile = "cpu" * filename;
gpufile = "gpu" * filename;

nbc = size(f,2);

stropu = "";
strgpu = "";
for k = 1:nbc
    str1, str2 = gencodebou(filename * string(k), f[:,k], xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time);
    stropu = stropu * str1;
    strgpu = strgpu * str2;
end

tmp = "template <typename T> void " * opufile;
tmp = tmp * "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";
tmp = tmp * "{\n";
for k = 1:nbc
    if k == 1
        tmp = tmp * "\tif (ib == " * string(k) * ")\n";
    else
        tmp = tmp * "\telse if (ib == " * string(k) * ")\n";
    end
    tmp = tmp * "\t\t" * opufile * string(k) * "(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, ng, nc, ncu, nd, ncx, nco, ncw);\n";
end
tmp = tmp * "}\n\n";

tmp = tmp * "template void " * opufile;
tmp = tmp * "(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int);\n";
tmp = tmp * "template void " * opufile;
tmp = tmp * "(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int);\n";
stropu = stropu * tmp;

ioopu = open(foldername * "/" * opufile * ".cpp", "w");
write(ioopu, stropu);
close(ioopu);

tmp = "template <typename T> void " * gpufile;
tmp = tmp * "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";
tmp = tmp * "{\n";
tmp = tmp * "\tint blockDim = 256;\n";
tmp = tmp * "\tint gridDim = (ng + blockDim - 1) / blockDim;\n";
tmp = tmp * "\tgridDim = (gridDim>1024)? 1024 : gridDim;\n";
for k = 1:nbc
    if k == 1
        tmp = tmp * "\tif (ib == " * string(k) * ")\n";
    elseif k<nbc
        tmp = tmp * "\telse if (ib == " * string(k) * ")\n";
    else
        tmp = tmp * "\telse if (ib == " * string(k) * ")\n";
    end
    tmp = tmp * "\t\t"  * gpufile * string(k) * "<<<gridDim, blockDim>>>(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, ng, nc, ncu, nd, ncx, nco, ncw);\n";
end
tmp = tmp * "}\n\n";

tmp = tmp * "template void " * gpufile;
tmp = tmp * "(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int);\n";
tmp = tmp * "template void " * gpufile;
tmp = tmp * "(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int);\n";
strgpu = strgpu * tmp;

iogpu = open(foldername * "/" * gpufile * ".cu", "w");
write(iogpu, strgpu);
close(iogpu);

strcpu = replace(stropu, "opu" => "cpu");
strcpu = replace(strcpu, "for (int i = 0; i <ng; i++) {" => "#pragma omp parallel for\n\tfor (int i = 0; i <ng; i++) {");

iocpu = open(foldername * "/" * cpufile * ".cpp", "w");
write(iocpu, strcpu);
close(iocpu);

end
