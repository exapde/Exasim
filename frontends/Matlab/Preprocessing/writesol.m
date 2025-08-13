function writesol(app, master, mesh, dmd, filename, i)    

endian = 'native';

disp(['Writing initial solution into file ...']);
if app.mpiprocs>0
    fileID = fopen(filename + "sol" + string(i) + ".bin",'w');
else
    fileID = fopen(filename + "sol" + ".bin",'w');
end

ndims = zeros(12,1);           
ndims(1) = length(dmd.elempart); % number of elements
ndims(2) = sum(dmd.facepartpts); % number of faces
ndims(3) = mesh.nfe; 
ndims(4) = master.npe;
ndims(5) = master.npf;            
ndims(6) = app.nc;
ndims(7) = app.ncu;
ndims(8) = app.ncq;
ndims(9) = app.ncw;
ndims(10) = app.nco;
ndims(11) = app.nch;
ndims(12) = app.ncx;

nsize = zeros(20,1);
nsize(1) = length(ndims(:));
if isfield(mesh, 'xdg')        
    nsize(2) = numel(mesh.xdg(:,:,dmd.elempart));
end
if isfield(mesh, 'udg')        
    nsize(3) = numel(mesh.udg(:,:,dmd.elempart));
end
if isfield(mesh, 'vdg')        
    nsize(4) = numel(mesh.vdg(:,:,dmd.elempart));
end
if isfield(mesh, 'wdg')        
    nsize(5) = numel(mesh.wdg(:,:,dmd.elempart));
end
if isfield(mesh, 'uhat')   
    nsize(6) = numel(mesh.uhat);
end

fwrite(fileID,length(nsize(:)),'double',endian);
fwrite(fileID,nsize(:),'double',endian);
fwrite(fileID,ndims(:),'double',endian);
if isfield(mesh, 'xdg')        
    fwrite(fileID,mesh.xdg(:,:,dmd.elempart),'double',endian);                
end
if isfield(mesh, 'udg')        
    fwrite(fileID,mesh.udg(:,:,dmd.elempart),'double',endian);                
end
if isfield(mesh, 'vdg')        
    fwrite(fileID,mesh.vdg(:,:,dmd.elempart),'double',endian);                
end
if isfield(mesh, 'wdg')        
    fwrite(fileID,mesh.wdg(:,:,dmd.elempart),'double',endian);                
end
if isfield(mesh, 'uhat')   
    fwrite(fileID,mesh.uhat,'double',endian);            
end

fclose(fileID);         
