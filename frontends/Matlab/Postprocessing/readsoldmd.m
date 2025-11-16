function [k, sol] = readsoldmd(base, dmd, nsteps, stepoffsets)

if nargin < 3
    nsteps = 1;
end
if nargin < 4
    stepoffsets = 0;
end

nproc = length(dmd);
k = zeros(nproc, 4);
for i = 1:nproc
    %fname = rankfile(base, i-1);
    fname = base + "_np" + string(i-1) + ".bin";
    
    fid = fopen(fname, 'r');
    if fid < 0, error('Cannot open file: %s', fname); end
    
    hdr = fread(fid, 3, 'double').';
    if numel(hdr) ~= 3
        fclose(fid); error('Header read failed for %s', fname);
    end
    n1 = hdr(1); n2 = hdr(2); n3 = hdr(3);
    N  = n1 * n2 * n3;
    
    info = dir(fname); 
    L = info.bytes/8;
    timesteps = (L-3)/N;

    k(i, :) = [n1, n2, n3, timesteps];

    fclose(fid);
end

k = round(k);

if sum(abs(k(:,1) - k(1,1))) ~= 0
    error("First dimension in binary files do not match.");
end
if sum(abs(k(:,2) - k(1,2))) ~= 0
    error("Second dimension in binary files do not match.");
end
if sum(abs(k(:,4) - k(1,4))) ~= 0
    error("Number of timesteps in binary files do not match.");
end

if nargout > 1 
    ne = sum(k(:,3));
    sol = zeros(n1,n2,ne,nsteps);    

    for n = 1:nproc
        %fname = rankfile(base, n-1);
        fname = base + "_np" + string(n-1) + ".bin";
        
        fid = fopen(fname, 'r');
        if fid < 0, error('Cannot open file: %s', fname); end

        hdr = fread(fid, 3, 'double').';
        if numel(hdr) ~= 3
            fclose(fid); error('Header read failed for %s', fname);
        end
        n1 = hdr(1); n2 = hdr(2); n3 = hdr(3);
        N  = n1 * n2 * n3;
        
        if (stepoffsets > 0) 
            fseek(fid, stepoffsets*N*8, 'cof');
        end
        
        elempart = dmd{n}.elempart(1:n3);        
        for i = 1:nsteps
            tm = fread(fid, N, 'double');
            sol(:,:,elempart,i) = reshape(tm, [n1, n2, n3]);
        end

        fclose(fid);
    end

end

end

