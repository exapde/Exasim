function updatesolfiles(app,dgnodes,UDG,ODG,WDG,elempart)

endian = 'native';
mpiprocs = length(elempart);
for i = 1:mpiprocs            
    xdg = dgnodes(:,:,elempart{i});                        
    if isempty(UDG)==0
        udg = UDG(:,:,elempart{i});
    else
        udg = [];
    end
    if isempty(ODG)==0
        odg = ODG(:,:,elempart{i});
    else
        odg = [];
    end
    if isempty(WDG)==0
        wdg = WDG(:,:,elempart{i});
    else
        wdg = [];
    end

    disp(['Writing initial solution into file ' num2str(i) '...']); disp(' ');
    fileID1 = fopen([app.filename,'sol' num2str(i) '.bin'],'w');
    ndims = zeros(12,1);           
    %ndims(1) = length(elempart{i}); % number of elements
    %ndims(2) = length(dmd{i}.facepart); % number of faces
    %ndims(3) = size(t2f,2); % number of faces per element          
    %ndims(4) = size(UDG,1);
    %ndims(5) = master.npf(1);            
    ndims(6) = app.nc;
    ndims(7) = app.ncu;
    ndims(8) = app.ncq;
    ndims(9) = app.ncp;
    ndims(10) = app.nco;
    ndims(11) = app.nch;
    ndims(12) = app.ncx;

    nsize = zeros(20,1);
    nsize(1) = length(ndims(:));
    nsize(2) = length(xdg(:));
    nsize(3) = length(udg(:)); 
    nsize(4) = length(odg(:));    
    nsize(5) = length(wdg(:));    

    fwrite(fileID1,length(nsize(:)),'double',endian);
    fwrite(fileID1,nsize(:),'double',endian);
    fwrite(fileID1,ndims(:),'double',endian);
    fwrite(fileID1,xdg(:),'double',endian);
    fwrite(fileID1,udg(:),'double',endian);        
    fwrite(fileID1,odg(:),'double',endian);
    fwrite(fileID1,wdg(:),'double',endian);
    fclose(fileID1);         
end
