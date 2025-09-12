function stropu = gencodeelemface(filename, npm, ncase, foldername)

%foldername = "app";
opufile = "opu" + filename;
cpufile = "cpu" + filename;
gpufile = "gpu" + filename;

tmp = "";
for k = 1:npm
    tmp = tmp + "#include " + """" + opufile + string(k) +  ".cpp" + """" + "\n";
end
tmp = tmp + "\n";

if ncase==1
    tmp = tmp + "template <typename T> void " + opufile;
    tmp = tmp + "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";
    tmp = tmp + "{\n";
    for k = 1:npm
        if k == 1
            tmp = tmp + "\tif (modelnumber == " + string(k) + ")\n";    
        else            
            tmp = tmp + "\telse if (modelnumber == " + string(k) + ")\n";    
        end 
        tmp = tmp + "\t\t" + opufile + string(k) + "(f, xdg, udg, odg, wdg, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);\n";
    end
    tmp = tmp + "}\n\n";
    tmp = tmp + "template void " + opufile;
    tmp = tmp + "(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int);\n";
    tmp = tmp + "template void " + opufile;
    tmp = tmp + "(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int);\n";
elseif ncase==2
    tmp = tmp + "template <typename T> void " + opufile;
    tmp = tmp + "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne)\n";
    tmp = tmp + "{\n";
    for k = 1:npm
        if k == 1
            tmp = tmp + "\tif (modelnumber == " + string(k) + ")\n";    
        else            
            tmp = tmp + "\telse if (modelnumber == " + string(k) + ")\n";    
        end 
        tmp = tmp + "\t\t" + opufile + string(k) + "(f, xdg, udg, odg, wdg, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne);\n";
    end
    tmp = tmp + "}\n\n";
    tmp = tmp + "template void " + opufile;
    tmp = tmp + "(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int, int, int);\n";
    tmp = tmp + "template void " + opufile;
    tmp = tmp + "(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int, int, int);\n";   
elseif ncase==3
    tmp = tmp + "template <typename T> void " + opufile;
    tmp = tmp + "(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";
    tmp = tmp + "{\n";
    for k = 1:npm
        if k == 1
            tmp = tmp + "\tif (modelnumber == " + string(k) + ")\n";    
        else            
            tmp = tmp + "\telse if (modelnumber == " + string(k) + ")\n";    
        end 
        tmp = tmp + "\t\t" + opufile + string(k) + "(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ib, ng, nc, ncu, nd, ncx, nco, ncw);\n";
    end
    tmp = tmp + "}\n\n";
    tmp = tmp + "template void " + opufile;
    tmp = tmp + "(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int);\n";
    tmp = tmp + "template void " + opufile;
    tmp = tmp + "(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int);\n";
elseif ncase==4
    tmp = tmp + "template <typename T> void " + opufile;
    tmp = tmp + "(T *f, T *xdg, T *udg1, T *udg2,  T *odg1, T *odg2,  T *wdg1, T *wdg2,  T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)\n";
    tmp = tmp + "{\n";
    for k = 1:npm
        if k == 1
            tmp = tmp + "\tif (modelnumber == " + string(k) + ")\n";    
        else            
            tmp = tmp + "\telse if (modelnumber == " + string(k) + ")\n";    
        end 
        tmp = tmp + "\t\t" + opufile + string(k) + "(f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);\n";
    end
    tmp = tmp + "}\n\n";
    tmp = tmp + "template void " + opufile;
    tmp = tmp + "(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double time, int, int, int, int, int, int, int, int);\n";
    tmp = tmp + "template void " + opufile;
    tmp = tmp + "(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float time, int, int, int, int, int, int, int, int);\n";
elseif ncase==5    
    tmp = tmp + "template <typename T> void " + opufile;
    tmp = tmp + "(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)\n";
    tmp = tmp + "{\n";
    for k = 1:npm
        if k == 1
            tmp = tmp + "\tif (modelnumber == " + string(k) + ")\n";    
        else            
            tmp = tmp + "\telse if (modelnumber == " + string(k) + ")\n";    
        end 
        tmp = tmp + "\t\t" + opufile + string(k) + "(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);\n";
    end
    tmp = tmp + "}\n\n";
    tmp = tmp + "template void " + opufile;
    tmp = tmp + "(double *, double *, double *, double *, int, int, int, int, int, int);\n";
    tmp = tmp + "template void " + opufile;
    tmp = tmp + "(float *, float *, float *, float *, int, int, int, int, int, int);\n";
end

stropu = tmp;
fid = fopen(foldername + "/" + opufile + ".cpp", 'w');
fprintf(fid, char(stropu));
fclose(fid);

strcpu = strrep(stropu, 'opu', "cpu");
iocpu = fopen(foldername + "/" + cpufile + ".cpp", 'w');
fprintf(fid, char(strcpu));
fclose(iocpu);

strgpu = strrep(stropu, 'opu', "gpu");
strgpu = strrep(strgpu, 'cpp', "cu");
iogpu = fopen(foldername + "/" + gpufile + ".cu", 'w');
fprintf(fid, char(strgpu));
fclose(iogpu);

end

