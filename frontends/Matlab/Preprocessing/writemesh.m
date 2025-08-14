function writemesh(app, master, mesh, dmd, filename, i)    

endian = 'native';
if (app.mpiprocs>1) 
    fileID = fopen(filename + "mesh" + string(i) + ".bin",'w');
else
    fileID = fopen(filename + "mesh" + ".bin",'w');
end

ndims = zeros(20,1);
ndims(1) = size(mesh.p,1);
ndims(2) = length(dmd.elempart);
ndims(3) = sum(dmd.facepartpts);
ndims(4) = max(mesh.t(:));
ndims(5) = mesh.nfe;
ndims(6) = mesh.nbe;
ndims(7) = mesh.neb;
ndims(8) = mesh.nbf;
ndims(9) = mesh.nfb;

nsize = zeros(50,1);
nsize(1) = length(ndims(:));
nsize(2) = length(dmd.facecon(:));  
nsize(3) = length(mesh.eblks(:)); 
nsize(4) = length(mesh.fblks(:)); 
nsize(5) = length(dmd.nbsd(:)); 
nsize(6) = length(dmd.elemsend(:)); 
nsize(7) = length(dmd.elemrecv(:)); 
nsize(8) = length(dmd.elemsendpts(:)); 
nsize(9) = length(dmd.elemrecvpts(:));             
nsize(10) = length(dmd.elempart(:)); 
nsize(11) = length(dmd.elempartpts(:)); 
nsize(12) = length(mesh.cgelcon(:));  
nsize(13) = length(mesh.rowent2elem(:));  
nsize(14) = length(mesh.cgent2dgent(:));  
nsize(15) = length(mesh.colent2elem(:));                          
nsize(16) = length(mesh.rowe2f1(:));  
nsize(17) = length(mesh.cole2f1(:));  
nsize(18) = length(mesh.ent2ind1(:));                          
nsize(19) = length(mesh.rowe2f2(:));  
nsize(20) = length(mesh.cole2f2(:));  
nsize(21) = length(mesh.ent2ind2(:));  
if (app.hybrid > 0)
  nsize(22) = length(dmd.f2t(:));  
  nsize(23) = length(dmd.elemcon(:));  
  nsize(24) = length(master.perm(:));  
  nsize(25) = length(dmd.bf(:));     
  nsize(26) = length(app.cartgridpart);          
end

fwrite(fileID,length(nsize(:)),'double',endian);
fwrite(fileID,nsize(:),'double',endian);
fwrite(fileID,ndims(:),'double',endian);
fwrite(fileID,dmd.facecon(:),'double',endian);
fwrite(fileID,mesh.eblks(:),'double',endian);
fwrite(fileID,mesh.fblks(:),'double',endian);
fwrite(fileID,dmd.nbsd(:),'double',endian);
fwrite(fileID,dmd.elemsend(:),'double',endian);
fwrite(fileID,dmd.elemrecv(:),'double',endian);
fwrite(fileID,dmd.elemsendpts(:),'double',endian);
fwrite(fileID,dmd.elemrecvpts(:),'double',endian);
fwrite(fileID,dmd.elempart(:),'double',endian);
fwrite(fileID,dmd.elempartpts(:),'double',endian);
fwrite(fileID,mesh.cgelcon(:),'double',endian);
fwrite(fileID,mesh.rowent2elem(:),'double',endian);
fwrite(fileID,mesh.cgent2dgent(:),'double',endian);
fwrite(fileID,mesh.colent2elem(:),'double',endian);            
fwrite(fileID,mesh.rowe2f1(:),'double',endian);
fwrite(fileID,mesh.cole2f1(:),'double',endian);
fwrite(fileID,mesh.ent2ind1(:),'double',endian);                 
fwrite(fileID,mesh.rowe2f2(:),'double',endian);
fwrite(fileID,mesh.cole2f2(:),'double',endian);
fwrite(fileID,mesh.ent2ind2(:),'double',endian);        
if (app.hybrid > 0) 
  fwrite(fileID,dmd.f2t(:),'double',endian);        
  fwrite(fileID,dmd.elemcon(:),'double',endian);        
  fwrite(fileID,master.perm(:),'double',endian);        
  fwrite(fileID,dmd.bf(:),'double',endian);                
  fwrite(fileID,app.cartgridpart(:),'double',endian);               
end
fclose(fileID);            
