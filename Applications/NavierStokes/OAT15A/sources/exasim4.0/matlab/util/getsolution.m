function UDG = getsolution(filename,nproc,npe,nc,elempart)

if nproc==1
    fileID = fopen([filename '_np0.bin'],'r');
    UDG = fread(fileID,'double');
    fclose(fileID);
    UDG = reshape(UDG,npe,nc,[]);
else
    if nargin<5
        load dmd.mat;
    end
    nei = zeros(1,nproc);
    for i = 1:nproc
        nei(i) = length(elempart{i});
    end
    ne = sum(nei);
    
    for i = 1:nproc                
        fileID = fopen([filename '_np' num2str(i-1) '.bin'],'r');
        tm = fread(fileID,npe*nc*nei(i),'double');
        UDG(:,:,elempart{i}) = reshape(tm,[npe nc nei(i)]);                   
        fclose(fileID);
    end    
end


