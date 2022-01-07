function mesh = block_crs(mesh,hybrid)

if strcmp(hybrid,'hdg')
    ne    = size(mesh.t2f,1);
    nfe   = size(mesh.t2f,2); 
    il = zeros(nfe,nfe,ne);
    jl = zeros(nfe,nfe,ne);
    for i=1:ne    
        con = mesh.t2f(i,:);
        il(:,:,i) = repmat(con' , 1, nfe);
        jl(:,:,i) = repmat(con , nfe, 1);        
    end
    il(:,:,1);
    jl(:,:,1);
    sl = ones(nfe,nfe,ne);

elseif strcmp(hybrid,'edg')
    
    [nn,ne] = size(mesh.elcon);
    %ne = mesh.ne;
    il = zeros(nn,nn,ne);
    jl = zeros(nn,nn,ne);
    for i=1:ne
        con = mesh.elcon(:,i);
        il(:,:,i) = repmat(con ,1,nn);
        jl(:,:,i) = repmat(con',nn,1);        
    end    
    sl = ones(nn,nn,ne);
end

[rp, cj] = sparse_to_csr(sparse(il(:),jl(:),sl(:)));

ct = cj;
for i = 1:length(rp)-1
    k = rp(i)+1:rp(i+1);    
    ct(k) = [i; sort(setdiff(ct(k),i))]; % the first index is the block itself   
end

if length(ct) ~= length(cj)
    error('Matrix must have nonzero diagonal entries!');
else
    cj = ct;
end

tm = rp(2:end)-rp(1:end-1);
mesh.maxBlocksPerRow = max(tm);  
mesh.minBlocksPerRow = min(tm);  
mesh.cbsr_rowpts = rp;
mesh.cbsr_colind = cj;
mesh.cbsr_nrows = length(rp)-1;
mesh.cbsr_nblks = length(cj);
mesh.hybrid = hybrid;


















