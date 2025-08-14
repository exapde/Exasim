function UH = fetchuhat(dmd,pde,timestep)

nproc = length(dmd);
UH = cell(nproc,1);
for i = 1:nproc
    if nargin<3
      fileID = fopen(pde.buildpath + "/dataout/out_uhat_np" + string(i-1) + ".bin",'r');
    else
      fileID = fopen(pde.buildpath + "/dataout/out_uhat_t" + string(timestep) + "_np" + string(i-1) + ".bin",'r');
    end
    ncu = pde.ncu;
    nf = size(dmd{i}.f2t,2);
    npf = size(dmd{i}.elemcon,1);
    UH{i} = reshape(fread(fileID,'double'), [ncu npf nf]);    
end

