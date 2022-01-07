function writemaster(master,filename)

ndims = zeros(20,1);
ndims[1] = master.nd;
ndims[2] = master.elemtype;
ndims[3] = master.nodetype;
ndims[4] = master.porder;
ndims[5] = master.pgauss;
ndims[6] = master.npe;
ndims[7] = master.npf;
ndims[8] = master.nge;
ndims[9] = master.ngf;
ndims[10] = length(master.xp1d);
ndims[11] = length(master.gp1d);

nsize = zeros(22,1);
nsize[1] = length(ndims[:]);
nsize[2] = length(master.shapegt[:]);
nsize[3] = length(master.shapegw[:]);
nsize[4] = length(master.shapfgt[:]);
nsize[5] = length(master.shapfgw[:]);
nsize[6] = length(master.shapent[:]);
nsize[7] = length(master.shapen[:]);
nsize[8] = length(master.shapfnt[:]);
nsize[9] = length(master.shapfn[:]);
nsize[10] = length(master.xpe[:]);
nsize[11] = length(master.gpe[:]);
nsize[12] = length(master.gwe[:]);
nsize[13] = length(master.xpf[:]);
nsize[14] = length(master.gpf[:]);
nsize[15] = length(master.gwf[:]);
nsize[16] = length(master.shap1dgt[:]);
nsize[17] = length(master.shap1dgw[:]);
nsize[18] = length(master.shap1dnt[:]);
nsize[19] = length(master.shap1dn[:]);
nsize[20] = length(master.xp1d[:]);
nsize[21] = length(master.gp1d[:]);
nsize[22] = length(master.gw1d[:]);

print("Writing master into file...\n");
fileID = open(filename,"w");

# write master structure to files
write(fileID,Float64(length(nsize[:])));
write(fileID,Float64.(nsize[:]));
write(fileID,Float64.(ndims[:]));
write(fileID,master.shapegt[:]);
write(fileID,master.shapegw[:]);
write(fileID,master.shapfgt[:]);
write(fileID,master.shapfgw[:]);
write(fileID,master.shapent[:]);
write(fileID,master.shapen[:]);
write(fileID,master.shapfnt[:]);
write(fileID,master.shapfn[:]);
write(fileID,master.xpe[:]);
write(fileID,master.gpe[:]);
write(fileID,master.gwe[:]);
write(fileID,master.xpf[:]);
write(fileID,master.gpf[:]);
write(fileID,master.gwf[:]);
write(fileID,master.shap1dgt[:]);
write(fileID,master.shap1dgw[:]);
write(fileID,master.shap1dnt[:]);
write(fileID,master.shap1dn[:]);
write(fileID,master.xp1d[:]);
write(fileID,master.gp1d[:]);
write(fileID,master.gw1d[:]);

close(fileID);

end
