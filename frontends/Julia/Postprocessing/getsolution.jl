function getsolution(filename,dmd,npe)

nproc = length(dmd);
if nproc==1
    fn = filename * "_np0.bin";
    UDG = reinterpret(Float64,read(fn));
    ne = length(dmd[1].elempart[:]);
    nc = Int64(round(length(UDG)/(npe*ne)));
    UDG = reshape(UDG,(npe,nc,ne));
else
    nei = zeros(Int,nproc);
    for i = 1:nproc
        nei[i] = sum(dmd[i].elempartpts[1:2]);
    end
    ne = sum(nei);

    fn = filename * "_np0.bin";
    tmp = reinterpret(Float64,read(fn));
    nc = Int64(round(length(tmp)/(npe*nei[1])));

    UDG = zeros(npe,nc,ne);
    for i = 1:nproc
        elempart = dmd[i].elempart[1:nei[i]];
        fn = filename * "_np" * string(i-1) * ".bin"
        UDG[:,:,elempart]  = reshape(reinterpret(Float64,read(fn)),(npe,nc,nei[i]));
    end
end

return UDG

end
