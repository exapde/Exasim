function writemesh(mesh,filename,ifile)


nbe = size(mesh.eblks,2);
nbf = size(mesh.fblks,2);
neb = maximum(eblks[2,:] .- eblks[1,:])+1;
nfb = maximum(fblks[2,:] .- fblks[1,:])+1;
# nfe = size(mesh.t2f,1);
# nv = maximum(mesh.t[:]));

ndims = zeros(20,1);
ndims[1] = mesh.nd;
ndims[2] = mesh.ne;
ndims[3] = mesh.nf;
ndims[4] = mesh.nv;
ndims[5] = mesh.nfe;
ndims[6] = nbe;
ndims[7] = neb;
ndims[8] = nbf;
ndims[9] = nfb;

nsize = zeros(30,1);
nsize[1] = length(ndims);
nsize[2] = length(mesh.facecon);
nsize[3] = length(mesh.eblks);
nsize[4] = length(mesh.fblks);
nsize[5] = length(mesh.nbsd);
nsize[6] = length(mesh.elemsend);
nsize[7] = length(mesh.elemrecv);
nsize[8] = length(mesh.elemsendpts);
nsize[9] = length(mesh.elemrecvpts);
nsize[10] = length(mesh.elempart);
nsize[11] = length(mesh.elempartpts);
nsize[12] = length(mesh.cgelcon);
nsize[13] = length(mesh.rowent2elem);
nsize[14] = length(mesh.cgent2dgent);
nsize[15] = length(mesh.colent2elem);
nsize[16] = length(mesh.rowe2f1);
nsize[17] = length(mesh.cole2f1);
nsize[18] = length(mesh.ent2ind1);
nsize[19] = length(mesh.rowe2f2);
nsize[20] = length(mesh.cole2f2);
nsize[21] = length(mesh.ent2ind2);

if ifile>0
    display(["Writing mesh into file " * string(ifile)]);
    fileID = open([filename * string(ifile) * ".bin"],"w");
else
    display(["Writing mesh into file "]);
    fileID = open([filename * ".bin"],"w");
end

write(fileID,Float64(length(nsize[:])));
write(fileID,Float64.(nsize[:]));

if nsize[1]>0
    write(fileID,Float64.(ndims[:]));
end
#reshape(permute(mesh.facecon,[2 1 3]),[2*master.npf(1)*ndims[3) 1])
if nsize[2]>0
    write(fileID,Float64.(mesh.facecon[:]));
end
if nsize[3]>0
    write(fileID,Float64.(mesh.eblks[:]));
end
if nsize[4]>0
    write(fileID,Float64.(mesh.fblks[:]));
end
if nsize[5]>0
    write(fileID,Float64.(mesh.nbsd[:]));
end
if nsize[6]>0
    write(fileID,Float64.(mesh.elemsend[:]));
end
if nsize[7]>0
    write(fileID,Float64.(mesh.elemrecv[:]));
end
if nsize[8]>0
    write(fileID,Float64.(mesh.elemsendpts[:]));
end
if nsize[9]>0
    write(fileID,Float64.(mesh.elemrecvpts[:]));
end
if nsize[10]>0
    write(fileID,Float64.(mesh.elempart[:]));
end
if nsize[11]>0
    write(fileID,Float64.(mesh.elempartpts[:]));
end
if nsize[12]>0
    write(fileID,Float64.(mesh.cgelcon[:]-1));
end
if nsize[13]>0
    write(fileID,Float64.(mesh.rowent2elem[:]));
end
if nsize[14]>0
    write(fileID,Float64.(mesh.cgent2dgent[:]-1));
end
if nsize[15]>0
    write(fileID,Float64.(mesh.colent2elem[:]-1));
end
if nsize[16]>0
    write(fileID,Float64.(mesh.rowe2f1[:]));
end
if nsize[17]>0
    write(fileID,Float64.(mesh.cole2f1[:]-1));
end
if nsize[18]>0
    write(fileID,Float64.(mesh.ent2ind1[:]-1));
end
if nsize[19]>0
    write(fileID,Float64.(mesh.rowe2f2[:]));
end
if nsize[20]>0
    write(fileID,Float64.(mesh.cole2f2[:]-1));
end
if nsize[21]>0
    write(fileID,Float64.(mesh.ent2ind2[:]-1));
end

close(fileID);

end
