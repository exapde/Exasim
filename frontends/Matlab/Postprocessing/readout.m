function sol = readout(base, elempartpts, elempart, npe)
nprocs = length(elempartpts);
Ne = 0;
Ni = zeros(nprocs,1);
for i = 1:nprocs
    Nj = sum(elempartpts{i}(1:2));
    Ne = Ne + Nj;
    Ni(i) = Nj;
end
filesol = base + string(0) + ".bin";

udg = readbin(filesol);
N1 = Ni(1);
nc = numel(udg)/(npe*N1);
ncu = nc;
sol = zeros(npe, ncu, Ne);
for i = 1:nprocs
    Nj = Ni(i);
    filesol = base + string(i-1) + ".bin";
    udg = readbin(filesol);
    if numel(udg) == npe*nc*Nj
        tm = reshape(udg, npe, nc, Nj);
        sol(:,:,elempart{i}(1:Nj)) = tm(:,1:ncu,:);
    else
        disp([i npe nc Nj])
    end
end
end







