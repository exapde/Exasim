function strkk = hdggencodefext(filename, f, xdg, udg, odg, wdg, uhg, nlg, uext, tau, uinf, param, time, foldername)

strkk = "";
nbc = size(f,2);
for k = 1:nbc
    str1 = gencodefext(filename + string(k), f(:,k), xdg, udg, odg, wdg, uhg, nlg, uext, tau, uinf, param, time);
    strkk = strkk + str1;    
end

cpufile = "Hdg" + filename;
tmp = "(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* uext, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n";
tmp = "void " + cpufile + tmp;
tmp = tmp + "{\n";
for k = 1:nbc
    if k == 1
        tmp = tmp + "\tif (ib == " + string(k) + ")\n";    
    else            
        tmp = tmp + "\telse if (ib == " + string(k) + ")\n";    
    end 
    tmp = tmp + "\t\t" + cpufile + string(k) + "(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, uext, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);\n";
end
tmp = tmp + "}\n\n";

strkk = strkk + tmp;
fid = fopen(foldername + "/" + cpufile + ".cpp", 'w');
fprintf(fid, char(strkk));
fclose(fid);

end

function strkk = gencodefext(filename, f, xdg, udg, odg, wdg, uhg, nlg, uext, tau, uinf, param, time)

cpufile = "Hdg" + filename;
tmp = "(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* uext, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)\n";
str = "\tKokkos::parallel_for(" + '"' + filename  + '"' + ", ng, KOKKOS_LAMBDA(const size_t i) {\n";
strkk = "void " + cpufile;
strkk = strkk + tmp + "{\n";
strkk = strkk + str;

fstr = string(f(:));
str = "";
str = varsassign(str, "param", length(param), 0, fstr);
str = varsassign(str, "uinf", length(uinf), 0, fstr);
str = varsassign(str, "tau", length(tau), 0, fstr);
str = varsassign(str, "xdg", length(xdg), 1, fstr);
str = varsassign(str, "udg", length(udg), 1, fstr);
str = varsassign(str, "uhg", length(uhg), 1, fstr);
str = varsassign(str, "odg", length(odg), 1, fstr);
str = varsassign(str, "wdg", length(wdg), 1, fstr);
str = varsassign(str, "nlg", length(nlg), 1, fstr);
str = varsassign(str, "uext", length(uext), 1, fstr);
str = symsassign2(str, f, udg, wdg, uhg);

strkk = strkk + str + "\t});\n" + "}\n\n"; 
strkk = strrep(strkk, "T ", "dstype ");

end



