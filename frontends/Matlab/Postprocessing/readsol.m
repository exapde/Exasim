function [n1, n2, n3, timesteps, sol] = readsol(base, nsteps, stepoffsets)

if nargin < 2
    nsteps = 1;
end
if nargin < 3
    stepoffsets = 0;
end

fname = rankfile(base, 0);

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

if nargout > 4 
    if (stepoffsets > 0) 
        fseek(fid, stepoffsets*N*8, 'cof');
    end
    sol = zeros(n1,n2,n3,nsteps);    
    for i = 1:nsteps
        tm = fread(fid, N, 'double');
        sol(:,:,:,i) = reshape(tm, [n1, n2, n3]);
    end
end

fclose(fid);

end
