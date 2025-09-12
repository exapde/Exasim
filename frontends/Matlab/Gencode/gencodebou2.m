function stropu = gencodebou2(filename, npm, foldername)

%foldername = "app";
opufile = "opu" + filename;
cpufile = "cpu" + filename;
gpufile = "gpu" + filename;

tmp = "";
for k = 1:npm
    tmp = tmp + "#include " + """" + opufile + string(k) +  ".cpp" + """" + "\n";
end
tmp = tmp + "\n";

tmp = tmp + "template <typename T> void " + opufile;
tmp = tmp + "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumner, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";
tmp = tmp + "{\n";
for k = 1:npm
    if k == 1
        tmp = tmp + "\tif (modelnumner == " + string(k) + ")\n";    
    else            
        tmp = tmp + "\telse if (modelnumner == " + string(k) + ")\n";    
    end 
    tmp = tmp + "\t\t" + opufile + string(k) + "(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, ib, ng, nc, ncu, nd, ncx, nco, ncw);\n";
end
tmp = tmp + "}\n\n";

tmp = tmp + "template void " + opufile;
tmp = tmp + "(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int);\n";
tmp = tmp + "template void " + opufile;
tmp = tmp + "(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int);\n";
stropu = tmp;

fid = fopen(foldername + "/" + opufile + ".cpp", 'w');
fprintf(fid, char(stropu));
fclose(fid);

strcpu = strrep(stropu, 'opu', "cpu");
iocpu = fopen(foldername + "/" + cpufile + ".cpp", 'w');
fprintf(fid, char(strcpu));
fclose(iocpu);

strgpu = strrep(stropu, 'opu', "gpu");
iogpu = fopen(foldername + "/" + gpufile + ".cu", 'w');
fprintf(fid, char(strgpu));
fclose(iogpu);

end
