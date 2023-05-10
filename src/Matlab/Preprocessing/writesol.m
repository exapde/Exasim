function [] = writesol(app,mesh,master,dmd)
    if app.modelnumber==0
        strn = "";
    else
        strn = num2str(app.modelnumber);
    end
    
%     if  ~exist(char("datain" + strn), 'dir')
%         mkdir(char("datain" + strn));
%     end
%     if  ~exist(char("dataout" + strn), 'dir')
%         mkdir(char("dataout" + strn));
%     end
    filename = "datain" + strn + "/";
    endian = 'native';
   
    mpiprocs = app.mpiprocs;
    
    for i = 1:mpiprocs            
        disp(['Writing initial solution into file ' num2str(i) '...']);
        if mpiprocs>1
            fileID1 = fopen(filename + "sol" + string(i) + ".bin",'w');
        else
            fileID1 = fopen(filename + "sol" + ".bin",'w');
        end
        xdg = mesh.dgnodes(:,:,dmd{i}.elempart);
        ndims = zeros(12,1);           
        ndims(1) = length(dmd{i}.elempart); % number of elements
        ndims(2) = sum(dmd{i}.facepartpts); % number of faces
        ndims(3) = size(master.perm,2); % number of faces per element          
        ndims(4) = master.npe;
        ndims(5) = master.npf;            
        ndims(6) = app.nc;
        ndims(7) = app.ncu;
        ndims(8) = app.ncq;
        ndims(9) = app.ncw;
        ndims(10) = app.nco;
        ndims(11) = app.nch;
        ndims(12) = app.ncx;
        app.nce=10;
        ndims(13) = app.nce;
    %     disp(master.perm);
    
        nsize = zeros(20,1);
        nsize(1) = length(ndims(:));
        nsize(2) = length(xdg(:));
        %nsize(3) = length(udg(:)); 
    %     nsize(4) = length(odg(:));    
    %     nsize(5) = length(wdg(:));    
        if isfield(mesh, 'udg')        
            nsize(3) = numel(mesh.udg(:,:,dmd{i}.elempart));
        end
        if isfield(mesh, 'vdg')        
            nsize(4) = numel(mesh.vdg(:,:,dmd{i}.elempart));
        end
        if isfield(mesh, 'wdg')        
            nsize(5) = numel(mesh.wdg(:,:,dmd{i}.elempart));
        end
        if isfield(mesh, 'dudg')        
            nsize(6) = numel(mesh.dudg(:,:,dmd{i}.elempart));
        end
    %     disp(nsize)
    
        fwrite(fileID1,length(nsize(:)),'double',endian);
        fwrite(fileID1,nsize(:),'double',endian);
        fwrite(fileID1,ndims(:),'double',endian);
        fwrite(fileID1,xdg(:),'double',endian);
        %fwrite(fileID1,udg(:),'double',endian);        
    %     fwrite(fileID1,odg(:),'double',endian);
    %     fwrite(fileID1,wdg(:),'double',endian);
        if isfield(mesh, 'udg')        
            fwrite(fileID1,mesh.udg(:,:,dmd{i}.elempart),'double',endian);                
        end
        if isfield(mesh, 'vdg')        
            fwrite(fileID1,mesh.vdg(:,:,dmd{i}.elempart),'double',endian);                
        end
        if isfield(mesh, 'wdg')        
            fwrite(fileID1,mesh.wdg(:,:,dmd{i}.elempart),'double',endian);                
        end
        if isfield(mesh, 'dudg')        
            fwrite(fileID1,mesh.dudg(:,:,dmd{i}.elempart),'double',endian);    
        end
        fclose(fileID1);         
    
        % divide elements and faces into blocks
    end


