function UDG = getsolution(filename,dmd,npe)

nproc = length(dmd);
if nproc==1
    fileID = fopen([filename '_np0.bin'],'r');
    UDG = fread(fileID,'double');
    fclose(fileID); 
    ne = length(dmd{1}.elempart(:));
    nc = numel(UDG)/(npe*ne);
    UDG = reshape(UDG,npe,nc,ne);
else
    nei = zeros(1,nproc);
    for i = 1:nproc
        nei(i) = sum(dmd{i}.elempartpts(1:2));
    end
    ne = sum(nei);
        
    fileID = fopen([filename '_np' num2str(0) '.bin'],'r');        
    tmp = fread(fileID,'double');   
    nc = numel(tmp)/(npe*nei(1));
    fclose(fileID);
    
    UDG = zeros(npe,nc,ne);    
    for i = 1:nproc        
        elempart = dmd{i}.elempart(1:nei(i));
        fileID = fopen([filename '_np' num2str(i-1) '.bin'],'r');        
        UDG(:,:,elempart) = reshape(fread(fileID,'double'),[npe nc nei(i)]);   
        fclose(fileID);
    end    
end


