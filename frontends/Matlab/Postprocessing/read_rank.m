function [n1, n2, n3, nsteps, payload] = read_rank(fname)
fid = fopen(fname, 'r');
if fid < 0, error('Cannot open file: %s', fname); end

hdr = fread(fid, 3, 'double').';
if numel(hdr) ~= 3
    fclose(fid); error('Header read failed for %s', fname);
end
n1 = hdr(1); n2 = hdr(2); n3 = hdr(3);
N  = n1 * n2 * n3;

if nargout > 3 
    payload = fread(fid, 'double');
    fclose(fid);
    
    if mod(numel(payload), N) ~= 0
        error('Payload size is not a multiple of n1*n2*n3 in %s', fname);
    end
    nsteps = numel(payload) / N;
    
    % reshape to [n1 n2 n3 nsteps]
    payload = reshape(payload, [n1, n2, n3, nsteps]);
end

end
