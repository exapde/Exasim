function writeudg(base, sol, dmd)

nprocs = length(dmd);

for i = 1:nprocs
    filesol = base + string(i) + ".bin";
    udg = sol(:,:,dmd{i}.elempart);
    writeudg2file(filesol, udg);
end

end

function writeudg2file(filesol, udg)
    
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

tm(k1:k2) = udg(:);
writebin(filesol,tm);

end

