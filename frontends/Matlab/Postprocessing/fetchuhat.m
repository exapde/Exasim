function sol = fetchuhat(dmd,pde,timestep)
% GETSOLUTIONS  Assemble solution blocks from per-rank binary files.
%
% File format per rank:
%   header (3 doubles): [n1, n2, n3]
%   followed by data for all timesteps, each block of size n1*n2*n3 (double)

if nargin < 2
    nproc = 1;
else
    nproc = numel(dmd);
end

basename = pde.buildpath + "/dataout/outuhat";

if nproc == 1
    [n1, n2, n3, nsteps, block] = read_rank(rankfile(basename, 0));
    sol = block;                                     % [n1 n2 n3 nsteps]
else
    sol = cell(nproc,1);
    for r = 1:nproc
        [n1r, n2r, n3r, nsteps_r, sol{r}] = read_rank(rankfile(basename, r-1));        
    end
end

end

% function UH = fetchuhat(dmd,pde,timestep)
% 
% nproc = length(dmd);
% UH = cell(nproc,1);
% for i = 1:nproc
%     if nargin<3
%       fileID = fopen(pde.buildpath + "/dataout/out_uhat_np" + string(i-1) + ".bin",'r');
%     else
%       fileID = fopen(pde.buildpath + "/dataout/out_uhat_t" + string(timestep) + "_np" + string(i-1) + ".bin",'r');
%     end
%     ncu = pde.ncu;
%     nf = size(dmd{i}.f2t,2);
%     npf = size(dmd{i}.elemcon,1);
%     UH{i} = reshape(fread(fileID,'double'), [ncu npf nf]);    
% end
% 
