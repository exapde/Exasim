function xdg = readxdgmpi(base, nprocs, ne)

if nprocs == 1
    filesol = base + ".bin";
    xdg = readxdg(filesol);        
    return;
end

for i = 1:nprocs
    filesol = base + string(i) + ".bin";
    if i == 1
        xdg = readxdg(filesol);    
        xdg = xdg(:,:,ne(i));
    else
        xdgi = readxdg(filesol);    
        xdg = cat(3, xdg, xdgi(:,:,ne(i)));
    end
end

