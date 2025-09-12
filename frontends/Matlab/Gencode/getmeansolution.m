function UDG = getmeansolution(filename,dmd,npe)

nproc = length(dmd);
if nproc==1
    fileID = fopen([filename '_np0.bin'],'r');
    UDG = fread(fileID,'double');
    fclose(fileID);
    ne = length(dmd{1}.elempart(:));
    nc = (numel(UDG)-1)/(npe*ne);
    UDG = reshape(UDG(1:(npe*nc*ne))/UDG(end),npe,nc,ne);
else
    nei = zeros(1,nproc);
    for i = 1:nproc
        nei(i) = sum(dmd{i}.elempartpts(1:2));
    end
    ne = sum(nei);
        
    fileID = fopen([filename '_np' num2str(0) '.bin'],'r');        
    tmp = fread(fileID,'double');   
    nc = (numel(tmp)-1)/(npe*nei(1));
        
    UDG = zeros(npe,nc,ne);    
    for i = 1:nproc        
        elempart = dmd{i}.elempart(1:nei(i));
        fileID = fopen([filename '_np' num2str(i-1) '.bin'],'r');        
        tmp = fread(fileID,'double');         
        UDG(:,:,elempart) = reshape(tmp(1:end-1),[npe nc nei(i)])/tmp(end);           
        fclose(fileID);
    end    
end


