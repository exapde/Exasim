function writemaster(master,filename,endian)

disp('Writing master into file...'); 

fileID = fopen(filename,'w');

ndims = zeros(20,1);
ndims(1) = master.dim;
ndims(2) = master.elemtype;
ndims(3) = master.nodetype;
ndims(4) = master.porder;
ndims(5) = master.pgauss;
ndims(6) = master.npe;
ndims(7) = master.npf;
ndims(8) = master.nge;
ndims(9) = master.ngf;
ndims(10) = master.np1d;
ndims(11) = master.ng1d;

nsize = zeros(22,1);
nsize(1) = length(ndims(:));
nsize(2) = length(master.shapegt(:));  
nsize(3) = length(master.shapegw(:)); 
nsize(4) = length(master.shapfgt(:)); 
nsize(5) = length(master.shapfgw(:));
nsize(6) = length(master.shapent(:));  
nsize(7) = length(master.shapen(:)); 
nsize(8) = length(master.shapfnt(:)); 
nsize(9) = length(master.shapfn(:));
nsize(10) = length(master.xpe(:)); 
nsize(11) = length(master.gpe(:)); 
nsize(12) = length(master.gwe(:)); 
nsize(13) = length(master.xpf(:)); 
nsize(14) = length(master.gpf(:)); 
nsize(15) = length(master.gwf(:)); 
nsize(16) = length(master.shap1dgt(:)); 
nsize(17) = length(master.shap1dgw(:)); 
nsize(18) = length(master.shap1dnt(:)); 
nsize(19) = length(master.shap1dn(:));
nsize(20) = length(master.xp1d(:)); 
nsize(21) = length(master.gp1d(:)); 
nsize(22) = length(master.gw1d(:)); 

% write master structure to files
fwrite(fileID,length(nsize(:)),'double',endian);
fwrite(fileID,nsize(:),'double',endian);
fwrite(fileID,ndims(:),'double',endian);
fwrite(fileID,master.shapegt(:),'double',endian);
fwrite(fileID,master.shapegw(:),'double',endian);
fwrite(fileID,master.shapfgt(:),'double',endian);
fwrite(fileID,master.shapfgw(:),'double',endian);
fwrite(fileID,master.shapent(:),'double',endian);
fwrite(fileID,master.shapen(:),'double',endian);
fwrite(fileID,master.shapfnt(:),'double',endian);
fwrite(fileID,master.shapfn(:),'double',endian);
fwrite(fileID,master.xpe(:),'double',endian);
fwrite(fileID,master.gpe(:),'double',endian);
fwrite(fileID,master.gwe(:),'double',endian);
fwrite(fileID,master.xpf(:),'double',endian);
fwrite(fileID,master.gpf(:),'double',endian);
fwrite(fileID,master.gwf(:),'double',endian);
fwrite(fileID,master.shap1dgt(:),'double',endian);
fwrite(fileID,master.shap1dgw(:),'double',endian);
fwrite(fileID,master.shap1dnt(:),'double',endian);
fwrite(fileID,master.shap1dn(:),'double',endian);
fwrite(fileID,master.xp1d(:),'double',endian);
fwrite(fileID,master.gp1d(:),'double',endian);
fwrite(fileID,master.gw1d(:),'double',endian);
fclose(fileID);
