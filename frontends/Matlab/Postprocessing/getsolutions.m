function sol = getsolutions(basename, dmd)
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

if nproc == 1
    [n1, n2, n3, nsteps, block] = read_rank(rankfile(basename, 0));
    sol = block;                                     % [n1 n2 n3 nsteps]
else
    % how many elements each rank contributes (like your original)
    nei = zeros(1, nproc);
    for r = 1:nproc
        nei(r) = sum(dmd{r}.elempartpts(1:2));
    end
    net = sum(nei);

    % read rank 0 to size global arrays
    [n1, n2, ~, nsteps, block0] = read_rank(rankfile(basename, 0));
    sol = zeros(n1, n2, net, nsteps);

    % place rank 0
    elempart = dmd{1}.elempart(1:nei(1));
    sol(:,:,elempart,:) = block0;

    % place ranks 1..nproc-1
    k = nei(1);
    for r = 2:nproc
        [n1r, n2r, n3r, nsteps_r, blockr] = read_rank(rankfile(basename, r-1));
        % (Optional) light sanity checks; comment out if layouts can differ:
        % if n1r~=n1 || n2r~=n2 || nsteps_r~=nsteps
        %     error('Rank %d has incompatible sizes.', r-1);
        % end
        elempart = dmd{r}.elempart(1:nei(r));
        sol(:,:,elempart,:) = blockr;
        k = k + nei(r);
    end
end

end

% ----------------- helpers -----------------

function fname = rankfile(base, r)
% Build "basename_np<r>.bin"
fname = base + "_np" + string(r) + ".bin";
end

function [n1, n2, n3, nsteps, data4d] = read_rank(fname)
fid = fopen(fname, 'r');
if fid < 0, error('Cannot open file: %s', fname); end

hdr = fread(fid, 3, 'double').';
if numel(hdr) ~= 3
    fclose(fid); error('Header read failed for %s', fname);
end
n1 = hdr(1); n2 = hdr(2); n3 = hdr(3);
N  = n1 * n2 * n3;

payload = fread(fid, 'double');
fclose(fid);

% if mod(numel(payload), N) ~= 0
%     error('Payload size is not a multiple of n1*n2*n3 in %s', fname);
% end
nsteps = floor(numel(payload) / N);

% reshape to [n1 n2 n3 nsteps]
data4d = reshape(payload(1:N*nsteps), [n1, n2, n3, nsteps]);
end

% function sol = getsolutions(basename, dmd)
% 
% if nargin<2
%     nproc = 1;
% else
%     nproc = length(dmd);
% end
% 
% if nproc==1    
%     fileID = fopen(basename + "_np0.bin", 'r');    
%     tmp = fread(fileID,'double');
%     fclose(fileID);
% 
%     n1 = tmp(1);
%     n2 = tmp(2);
%     n3 = tmp(3);
%     N = n1*n2*n3;    
%     nsteps = (length(tmp)-3)/N;
%     sol = zeros(n1, n2, n3, nsteps);
% 
%     m = 4; 
%     for i = 1:nsteps
%         sol(:,:,:,i) = reshape(tmp(m:m+N-1), [n1, n2, n3]);        
%         m = m + N;
%     end    
% else
%     nei = zeros(1,nproc);
%     for i = 1:nproc
%         nei(i) = sum(dmd{i}.elempartpts(1:2));
%     end
%     net = sum(nei);
% 
%     fileID = fopen(basename + "_np0.bin", 'r');    
%     tmp = fread(fileID,'double');
%     fclose(fileID);
% 
%     n1 = tmp(1);
%     n2 = tmp(2);
%     n3 = tmp(3);
%     N = n1*n2*n3;    
%     nsteps = (length(tmp)-3)/N;
% 
%     sol = zeros(n1, n2, net, nsteps);
% 
%     elempart = dmd{1}.elempart(1:nei(1));
%     m = 4;
%     for i = 1:nsteps
%         sol(:,:,elempart,i) = reshape(tmp(m:m+N-1), [n1, n2, n3]);        
%         m = m + N;
%     end    
% 
%     for n = 2:nproc
%         fileID = fopen(basename + "_np" + num2str(n-1) + ".bin", 'r');    
%         tmp = fread(fileID,'double');
%         fclose(fileID);
%         n1 = tmp(1);
%         n2 = tmp(2);
%         n3 = tmp(3);
%         N = n1*n2*n3;    
% 
%         elempart = dmd{n}.elempart(1:nei(n));
%         m = 4;
%         for i = 1:nsteps
%             sol(:,:,elempart,i) = reshape(tmp(m:m+N-1), [n1, n2, n3]);            
%             m = m + N;
%         end            
%     end
% end


% function [UDG, WDG, UHAT] = getsolutions(pde, dmd)
% 
% if nargin<2
%     nproc = 1;
% else
%     nproc = length(dmd);
% end
% 
% if nproc==1    
%     fileID = fopen(pde.buildpath + "/dataout/outsol_np0.bin", 'r');    
%     tmp = fread(fileID,'double');
%     fclose(fileID);
% 
%     npe = tmp(1);
%     nc = tmp(2);
%     ne = tmp(3);
%     ncw = tmp(5);
%     ncu = tmp(7);
%     npf = tmp(8);
%     nf = tmp(9);
%     N = tmp(10);    
%     nsteps = (length(tmp)-10)/N;
%     UDG = zeros(npe, nc, ne, nsteps);
%     WDG = zeros(npe, ncw, ne, nsteps);
%     UHAT = zeros(ncu, npf, nf, nsteps);
% 
%     n1 = npe*nc*ne;
%     n2 = npe*ncw*ne;
%     n3 = ncu*npf*nf;
%     m = 11; % Starting index for reshaping tmp into UDG    
%     for i = 1:nsteps
%         UDG(:,:,:,i) = reshape(tmp(m:m+n1-1), [npe, nc, ne]);
%         WDG(:,:,:,i) = reshape(tmp(m+n1:m+n1+n2-1), [npe, ncw, ne]);
%         UHAT(:,:,:,i) = reshape(tmp(m+n1+n2:m+n1+n2+n3-1), [ncu, npf, nf]);
%         m = m + N;
%     end    
% else
%     nei = zeros(1,nproc);
%     for i = 1:nproc
%         nei(i) = sum(dmd{i}.elempartpts(1:2));
%     end
%     net = sum(nei);
% 
%     nfi = zeros(1, nproc);
%     for n = 2:nproc
%         fileID = fopen(pde.buildpath + "/dataout/outsol_np" + num2str(n-1) + ".bin", 'r');    
%         tmp = fread(fileID, 10, 'double');
%         fclose(fileID);
%         nfi(n) = tmp(9);
%     end
%     nft = sum(nfi);
% 
%     fileID = fopen(pde.buildpath + "/dataout/outsol_np0.bin", 'r');    
%     tmp = fread(fileID,'double');
%     fclose(fileID);
% 
%     npe = tmp(1);
%     nc = tmp(2);
%     ne = tmp(3);
%     ncw = tmp(5);
%     ncu = tmp(7);
%     npf = tmp(8);
%     nf = tmp(9);
%     N = tmp(10);    
%     nsteps = (length(tmp)-10)/N;
% 
%     UDG = zeros(npe, nc, net, nsteps);
%     WDG = zeros(npe, ncw, net, nsteps);
%     UHAT = zeros(ncu, npf, nft, nsteps);
% 
%     elempart = dmd{1}.elempart(1:nei(1));
%     n1 = npe*nc*ne;
%     n2 = npe*ncw*ne;    
%     n3 = ncu*npf*nf;
%     m = 11;
%     for i = 1:nsteps
%         UDG(:,:,elempart,i) = reshape(tmp(m:m+n1-1), [npe, nc, ne]);
%         WDG(:,:,elempart,i) = reshape(tmp(m+n1:m+n1+n2-1), [npe, ncw, ne]);
%         UHAT(:,:,1:nf,i) = reshape(tmp(m+n1+n2:m+n1+n2+n3-1), [ncu, npf, nf]);
%         m = m + N;
%     end    
% 
%     k = nf;
%     for n = 2:nproc
%         fileID = fopen(pde.buildpath + "/dataout/outsol_np" + num2str(n-1) + ".bin", 'r');    
%         tmp = fread(fileID,'double');
%         fclose(fileID);
%         ne = tmp(3);
%         nf = tmp(9);
%         N = tmp(10);
% 
%         elempart = dmd{n}.elempart(1:nei(n));
%         n1 = npe*nc*ne;
%         n2 = npe*ncw*ne;    
%         n3 = ncu*npf*nf;
%         m = 11;
%         for i = 1:nsteps
%             UDG(:,:,elempart,i) = reshape(tmp(m:m+n1-1), [npe, nc, ne]);
%             WDG(:,:,elempart,i) = reshape(tmp(m+n1:m+n1+n2-1), [npe, ncw, ne]);
%             UHAT(:,:,k+1:k+nf,i) = reshape(tmp(m+n1+n2:m+n1+n2+n3-1), [ncu, npf, nf]);
%             m = m + N;
%         end    
%         k = k + nf;
%     end
% end
% 
% 
