function gencodeface(filename, f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time)

foldername = "app";
opufile = "opu" + filename;
cpufile = "cpu" + filename;
gpufile = "gpu" + filename;

nbc = size(f,2);

stropu = "";
strgpu = "";
for k = 1:nbc
    [str1, str2] = gencodebou(filename + string(k), f(:,k), xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time);
    stropu = stropu + str1;
    strgpu = strgpu + str2;
end

tmp = "template <typename T> void " + opufile;
tmp = tmp + "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";
tmp = tmp + "{\n";
for k = 1:nbc
    if k == 1
        tmp = tmp + "\tif (ib == " + string(k) + ")\n";    
    else            
        tmp = tmp + "\telse if (ib == " + string(k) + ")\n";    
    end 
    tmp = tmp + "\t\t" + opufile + string(k) + "(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);\n";
end
tmp = tmp + "}\n\n";

tmp = tmp + "template void " + opufile;
tmp = tmp + "(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int);\n";
tmp = tmp + "template void " + opufile;
tmp = tmp + "(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int);\n";
stropu = stropu + tmp;

fid = fopen(foldername + "/" + opufile + ".cpp", 'w');
fprintf(fid, char(stropu));
fclose(fid);

tmp = "template <typename T> void " + gpufile;
tmp = tmp + "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";
tmp = tmp + "{\n";
tmp = tmp + "\tint blockDim = 256;\n";
tmp = tmp + "\tint gridDim = (ng + blockDim - 1) / blockDim;\n";
tmp = tmp + "\tgridDim = (gridDim>1024)? 1024 : gridDim;\n";
for k = 1:nbc
    if k == 1
        tmp = tmp + "\tif (ib == " + string(k) + ")\n";    
    elseif k<nbc
        tmp = tmp + "\telse if (ib == " + string(k) + ")\n";    
    else            
        tmp = tmp + "\telse if (ib == " + string(k) + ")\n";    
    end 
    tmp = tmp + "\t\tkernel"  + gpufile + string(k) + "<<<gridDim, blockDim>>>(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);\n";
end
tmp = tmp + "}\n\n";

tmp = tmp + "template void " + gpufile;
tmp = tmp + "(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int);\n";
tmp = tmp + "template void " + gpufile;
tmp = tmp + "(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int);\n";
strgpu = strgpu + tmp;

% Enzyme
tmp = "#ifdef _ENZYME\n";
for k = 1:nbc
    if gpufile == "gpuUbou"
        tmp = tmp + "template <typename T> __global__ void kernelGrad" + gpufile  + string(k) + "Enzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *wdg, T *dwdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";
    else
        tmp = tmp + "template <typename T> __global__ void kernelGrad" + gpufile + string(k) + "Enzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *dodg, T *wdg, T *dwdg, T *uhg, T *duhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";
    end
    tmp = tmp + "{\n";
    tmp = tmp + "\t__enzyme_fwddiff" + gpufile + "((void*)device"  + gpufile + string(k) + "<T>,\n";
    tmp = tmp + "\t\t\t  enzyme_dup, f, df,\n";
    tmp = tmp + "\t\t\t enzyme_const, xg,\n";
    tmp = tmp + "\t\t\t enzyme_dup, udg, dudg,\n";
    if gpufile == "gpuUbou"
        tmp = tmp + "\t\t\t enzyme_const, odg,\n";
    else
        tmp = tmp + "\t\t\t enzyme_dup, odg, dodg,\n";
    end
    
    tmp = tmp + "\t\t\t enzyme_dup, wdg, dwdg,\n";
    if gpufile == "gpuUbou"
        tmp = tmp + "\t\t\t enzyme_const, uhg,\n";
    else
        tmp = tmp + "\t\t\t enzyme_dup, uhg, duhg,\n";
    end
    tmp = tmp + "\t\t\t enzyme_const, nlg,\n";
    tmp = tmp + "\t\t\t enzyme_const, tau,\n";
    tmp = tmp + "\t\t\t enzyme_const, uinf,\n";
    tmp = tmp + "\t\t\t enzyme_const, param,\n";
    tmp = tmp + "\t\t\t enzyme_const, time,\n";
    tmp = tmp + "\t\t\t enzyme_const, modelnumber,\n";
    tmp = tmp + "\t\t\t enzyme_const, ng,\n";
    tmp = tmp + "\t\t\t enzyme_const, nc,\n";
    tmp = tmp + "\t\t\t enzyme_const, ncu,\n";
    tmp = tmp + "\t\t\t enzyme_const, nd,\n";
    tmp = tmp + "\t\t\t enzyme_const, ncx,\n";
    tmp = tmp + "\t\t\t enzyme_const, nco,\n";
    tmp = tmp + "\t\t\t enzyme_const, ncw);\n";
    tmp = tmp + "}\n\n";
end


tmp = tmp + "template <typename T> void " + gpufile + "Enzyme";
if gpufile == "gpuUbou"
    tmp = tmp + "(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *wdg, T *dwdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";
else
    tmp = tmp + "(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *dodg, T *wdg, T *dwdg, T *uhg, T *duhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";
end
tmp = tmp + "{\n";
tmp = tmp + "\tint blockDim = 256;\n";
tmp = tmp + "\tint gridDim = (ng + blockDim - 1) / blockDim;\n";
tmp = tmp + "\tgridDim = (gridDim>1024)? 1024 : gridDim;\n";
for k = 1:nbc
    if k == 1
        tmp = tmp + "\tif (ib == " + string(k) + ")\n";
    elseif k<nbc
        tmp = tmp + "\telse if (ib == " + string(k) + ")\n";
    else
        tmp = tmp + "\telse if (ib == " + string(k) + ")\n";
    end
    if gpufile == "gpuUbou"
        tmp = tmp + "\t\tkernelGrad"  + gpufile + string(k) + "Enzyme<<<gridDim, blockDim>>>(f, df, xg, udg, dudg, odg, wdg, dwdg, uhg, nlg, tau, uinf, param, time, modelnumber, ib, ng, nc, ncu, nd, ncx, nco, ncw);\n";
    else
        tmp = tmp + "\t\tkernelGrad"  + gpufile + string(k) + "Enzyme<<<gridDim, blockDim>>>(f, df, xg, udg, dudg, odg, dodg, wdg, dwdg, uhg, duhg, nlg, tau, uinf, param, time, modelnumber, ib, ng, nc, ncu, nd, ncx, nco, ncw);\n";
    end
end
tmp = tmp + "}\n\n";

tmp = tmp + "template void " + gpufile + "Enzyme";
if gpufile == "gpuUbou"
    tmp = tmp + "(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int);\n";
else
    tmp = tmp + "(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int);\n";
end
tmp = tmp + "#endif";
strgpu = strgpu + tmp;
% end enzyme

fid = fopen(foldername + "/" + gpufile + ".cu", 'w');
fprintf(fid, char(strgpu));
fclose(fid);

strcpu = strrep(stropu, 'opu', "cpu");
strcpu = strrep(strcpu, "for (int i = 0; i <ng; i++) {", "#pragma omp parallel for\n\tfor (int i = 0; i <ng; i++) {");

iocpu = fopen(foldername + "/" + cpufile + ".cpp", 'w');
fprintf(fid, char(strcpu));
fclose(iocpu);

end
