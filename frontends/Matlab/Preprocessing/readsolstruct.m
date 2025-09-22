function sol = readsolstruct(filesol)
    
tm = readbin(filesol);
sz = tm(1);
k1 = 2;
k2 = k1+(sz)-1;
sol.nsize = tm(k1:k2);

k1 = k2+1;
k2 = k1+sol.nsize(1)-1;
sol.ndims = tm(k1:k2);

k1 = k2+1;
k2 = k1+sol.nsize(2)-1;
sol.xdg = tm(k1:k2);

k1 = k2+1;
k2 = k1+sol.nsize(3)-1;
sol.udg = tm(k1:k2);

k1 = k2+1;
k2 = k1+sol.nsize(4)-1;
sol.vdg = tm(k1:k2);

k1 = k2+1;
k2 = k1+sol.nsize(5)-1;
sol.wdg = tm(k1:k2);

% 
%         fileID1 = fopen(filename + "sol" + ".bin",'w');
%     end
%     ndims = zeros(12,1);           
%     ndims(1) = length(dmd{i}.elempart); % number of elements
%     ndims(2) = sum(dmd{i}.facepartpts); % number of faces
%     ndims(3) = size(master.perm,2); % number of faces per element          
%     ndims(4) = master.npe;
%     ndims(5) = master.npf;            
%     ndims(6) = app.nc;
%     ndims(7) = app.ncu;
%     ndims(8) = app.ncq;
%     ndims(9) = app.ncw;
%     ndims(10) = app.nco;
%     ndims(11) = app.nch;
%     ndims(12) = app.ncx;
% 
%     nsize = zeros(20,1);
%     nsize(1) = length(ndims(:));
%     nsize(2) = length(xdg(:));
%     %nsize(3) = length(udg(:)); 
% %     nsize(4) = length(odg(:));    
% %     nsize(5) = length(wdg(:));    
%     if isfield(mesh, 'udg')        
%         nsize(3) = numel(mesh.udg(:,:,dmd{i}.elempart));
%     end
%     if isfield(mesh, 'vdg')        
%         nsize(4) = numel(mesh.vdg(:,:,dmd{i}.elempart));
%     end
%     if isfield(mesh, 'wdg')        
%         nsize(5) = numel(mesh.wdg(:,:,dmd{i}.elempart));
%     end
% 
%     fwrite(fileID1,length(nsize(:)),'double',endian);
%     fwrite(fileID1,nsize(:),'double',endian);
%     fwrite(fileID1,ndims(:),'double',endian);
%     fwrite(fileID1,xdg(:),'double',endian);
%     if isfield(mesh, 'udg')        
%         fwrite(fileID1,mesh.udg(:,:,dmd{i}.elempart),'double',endian);                
%     end
%     if isfield(mesh, 'vdg')        
%         fwrite(fileID1,mesh.vdg(:,:,dmd{i}.elempart),'double',endian);                
%     end
%     if isfield(mesh, 'wdg')        
%         fwrite(fileID1,mesh.wdg(:,:,dmd{i}.elempart),'double',endian);                
%     end
%     fclose(fileID1);         

