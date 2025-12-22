function udgf = getudgf(base, nprocs, steps, stepoffsets)

udgf = [];
for i = 1:nprocs
    fname = base + "_np" + string(i-1) + ".bin";
    fid = fopen(fname, 'r');
    if fid < 0, error('Cannot open file: %s', fname); end
    hdr = fread(fid, 3, 'double').';
    if numel(hdr) ~= 3
        fclose(fid); error('Header read failed for %s', fname);
    end
    n1 = hdr(1); n2 = hdr(2); n3 = hdr(3);
    N = n1*n2*n3;
    if N > 0
        if (stepoffsets > 0)
          fseek(fid, stepoffsets*N*8, 'cof');
        end 
        tm = reshape(fread(fid, N, 'double'), [n1 n2 n3 steps]);
        if isempty(udgf)
            udgf = tm;
        else          
            udgf = cat(2, udgf, tm);
        end
    end
    fclose(fid);
end



