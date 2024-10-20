function writesol(filename,ifile,xdg, udg, odg, wdg)

if size(xdg,1)>0
    ncx = size(xdg,2);
else
    ncx = 0;
end
if size(udg,1)>0
    ncu = size(udg,2);
else
    ncu = 0;
end
if size(odg,1)>0
    nco = size(odg,2);
else
    nco = 0;
end
if size(wdg,1)>0
    ncw = size(wdg,2);
else
    ncw = 0;
end

ndims = zeros(12,1);
ndims[1] = ncx;
ndims[2] = ncu;
ndims[3] = nco;
ndims[4] = ncw;

nsize = zeros(20,1);
nsize[1] = length(ndims[:]);
nsize[2] = length(xdg[:]);
nsize[3] = length(udg[:]);
nsize[4] = length(odg[:]);
nsize[5] = length(wdg[:]);

if (ifile>0)
    print("Writing solution into file " * string(ifile) * "\n");
    fileID = open(string(filename, string(string(ifile), ".bin")),"w");
else
    print("Writing solution into file\n");
    fileID = open(string(filename, ".bin"),"w");
end

write(fileID,Float64(length(nsize[:])));
write(fileID,Float64.(nsize[:]));

if nsize[1]>0
    write(fileID,Float64.(ndims[:]));
end
if nsize[2]>0
    write(fileID,xdg[:]);
end
if nsize[3]>0
    write(fileID,udg[:]);
end
if nsize[4]>0
    write(fileID,odg[:]);
end
if nsize[5]>0
    write(fileID,wdg[:]);
end

close(fileID);

end
