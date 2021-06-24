function writesol2files(app,UDG,ODG,WDG,dgnodes,elempart,endian)

if nargin < 7; endianType = 0; end
if endianType == 0; endian = 'native';
elseif endianType == 1; endian = 'ieee-le';
elseif endianType == 2; endian = 'ieee-be';
end

if isfield(app,'filedir') == 0
    filename = app.appname;
else
    filename = app.filedir;
end

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
    fileID1 = fopen([filename,'sol' num2str(i) '.bin'],'w');
    ndims = zeros(12,1);           
    ndims(1) = length(elempart{i}); % number of elements
    ndims(2) = 1; % number of faces
    ndims(3) = 6; % number of faces per element          
    ndims(4) = size(dgnodes,1);
    ndims(5) = 9;            
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


