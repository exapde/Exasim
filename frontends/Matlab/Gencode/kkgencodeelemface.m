function stropu = kkgencodeelemface(filename, npm, ncase, foldername)

tmp = "";
for k = 1:npm
    tmp = tmp + "#include " + """" + filename + string(k) +  ".cpp" + """" + "\n";
end
tmp = tmp + "\n";

if ncase==1
    tmp = tmp + "void " + filename;
    tmp = tmp + "(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n";
    tmp = tmp + "{\n";
    for k = 1:npm
        if k == 1
            tmp = tmp + "\tif (modelnumber == " + string(k) + ")\n";    
        else            
            tmp = tmp + "\telse if (modelnumber == " + string(k) + ")\n";    
        end 
        tmp = tmp + "\t\t" + filename + string(k) + "(f, xdg, udg, odg, wdg, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);\n";
    end
    tmp = tmp + "}\n\n";
elseif ncase==2
    tmp = tmp + "void " + filename;
    tmp = tmp + "(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw, const int nce, const int npe, const int ne)\n";
    tmp = tmp + "{\n";
    for k = 1:npm
        if k == 1
            tmp = tmp + "\tif (modelnumber == " + string(k) + ")\n";    
        else            
            tmp = tmp + "\telse if (modelnumber == " + string(k) + ")\n";    
        end 
        tmp = tmp + "\t\t" + filename + string(k) + "(f, xdg, udg, odg, wdg, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne);\n";
    end
    tmp = tmp + "}\n\n";    
elseif ncase==3
    tmp = tmp + "void " + filename;
    tmp = tmp + "(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n";
    tmp = tmp + "{\n";
    for k = 1:npm
        if k == 1
            tmp = tmp + "\tif (modelnumber == " + string(k) + ")\n";    
        else            
            tmp = tmp + "\telse if (modelnumber == " + string(k) + ")\n";    
        end 
        tmp = tmp + "\t\t" + filename + string(k) + "(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ib, ng, nc, ncu, nd, ncx, nco, ncw);\n";
    end
    tmp = tmp + "}\n\n";
elseif ncase==4
    tmp = tmp + "void " + filename;
    tmp = tmp + "(dstype* f, const dstype* xdg, const dstype* udg1, const dstype* udg2,  const dstype* odg1, const dstype* odg2,  const dstype* wdg1, const dstype* wdg2,  const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n";
    tmp = tmp + "{\n";
    for k = 1:npm
        if k == 1
            tmp = tmp + "\tif (modelnumber == " + string(k) + ")\n";    
        else            
            tmp = tmp + "\telse if (modelnumber == " + string(k) + ")\n";    
        end 
        tmp = tmp + "\t\t" + filename + string(k) + "(f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);\n";
    end
    tmp = tmp + "}\n\n";
elseif ncase==5    
    tmp = tmp + "void " + filename;
    tmp = tmp + "(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)\n";
    tmp = tmp + "{\n";
    for k = 1:npm
        if k == 1
            tmp = tmp + "\tif (modelnumber == " + string(k) + ")\n";    
        else            
            tmp = tmp + "\telse if (modelnumber == " + string(k) + ")\n";    
        end 
        tmp = tmp + "\t\t" + filename + string(k) + "(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);\n";
    end
    tmp = tmp + "}\n\n";
elseif ncase==6    
    tmp = tmp + "void " + filename;
    tmp = tmp + "(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n";
    tmp = tmp + "{\n";
    for k = 1:npm
        if k == 1
            tmp = tmp + "\tif (modelnumber == " + string(k) + ")\n";    
        else            
            tmp = tmp + "\telse if (modelnumber == " + string(k) + ")\n";    
        end 
        tmp = tmp + "\t\t" + filename + string(k) + "(f, f_udg, f_wdg, xdg, udg, odg, wdg, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);\n";
    end
    tmp = tmp + "}\n\n";    
elseif ncase==7    
    tmp = tmp + "void " + filename;
    tmp = tmp + "(dstype* f, dstype* f_wdg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n";
    tmp = tmp + "{\n";
    for k = 1:npm
        if k == 1
            tmp = tmp + "\tif (modelnumber == " + string(k) + ")\n";    
        else            
            tmp = tmp + "\telse if (modelnumber == " + string(k) + ")\n";    
        end 
        tmp = tmp + "\t\t" + filename + string(k) + "(f, f_wdg, xdg, udg, odg, wdg, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);\n";
    end
    tmp = tmp + "}\n\n";     
elseif ncase==8    
    tmp = tmp + "void " + filename;
    tmp = tmp + "(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n";
    tmp = tmp + "{\n";
    for k = 1:npm
        if k == 1
            tmp = tmp + "\tif (modelnumber == " + string(k) + ")\n";    
        else            
            tmp = tmp + "\telse if (modelnumber == " + string(k) + ")\n";    
        end 
        tmp = tmp + "\t\t" + filename + string(k) + "(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ib, ng, nc, ncu, nd, ncx, nco, ncw);\n";
    end
    tmp = tmp + "}\n\n";  
elseif ncase==9    
    tmp = tmp + "void " + filename;
    tmp = tmp + "(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n";
    tmp = tmp + "{\n";
    for k = 1:npm
        if k == 1
            tmp = tmp + "\tif (modelnumber == " + string(k) + ")\n";    
        else            
            tmp = tmp + "\telse if (modelnumber == " + string(k) + ")\n";    
        end 
        tmp = tmp + "\t\t" + filename + string(k) + "(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ib, ng, nc, ncu, nd, ncx, nco, ncw);\n";
    end
    tmp = tmp + "}\n\n";                
end

stropu = tmp;
fid = fopen(foldername + "/" + filename + ".cpp", 'w');
fprintf(fid, char(stropu));
fclose(fid);

end

