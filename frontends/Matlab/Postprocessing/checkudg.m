function checkudg(base, sol, dmd)

nprocs = length(dmd);

for i = 1:nprocs
    filesol = base + string(i) + ".bin";
    udg = sol(:,:,dmd{i}.elempart);
    checkudg2file(filesol, udg);
end

end

function checkudg2file(filesol, udg)
    
tm = readbin(filesol);
sz = tm(1);
k1 = 2;
k2 = k1+(sz)-1;
nsize = tm(k1:k2);

k1 = k2+1;
k2 = k1+nsize(1)-1;

k1 = k2+1;
k2 = k1+nsize(2)-1;

k1 = k2+1;
k2 = k1+nsize(3)-1;

e = tm(k1:k2) - udg(:);
if max(abs(e(:))) > 1e-10
    error(filesol + " does not match udg!");
end

end

