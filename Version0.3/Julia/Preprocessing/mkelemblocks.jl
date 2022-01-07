function mkelemblocks(ne,ns::Int=512)

ns = min(ns,ne);
nb = Int(floor(ne/ns));  # number of blocks
na = Int(round(ne/nb)); # number of elements per block
nk = collect(1:na:ne);
nk = reshape(nk,1,length(nk));
if length(nk)==1
    nm = [1; ne];
else
    tm = 0*nk;
    tm[1:end-1] = nk[2:end].-1;
    tm[end] = ne;
    nm = [nk; tm];
    if (nm[2,end]-nm[1,end]) < (na/2)
        nm[2,end-1] = nm[2,end];
        nm = nm[:,1:(end-1)];
    end
end

nb = size(nm,2);
while minimum(nm[2,:]-nm[1,:])<(ns/2-1) || maximum(nm[2,:]-nm[1,:])>ns
    nb = nb+1;
    na = round(ne/nb); # number of elements per block
    nk = collect(1:na:ne);
    nk = reshape(nk,1,length(nk));
    tm = 0*nk;
    tm[1:end-1] = nk[2:end].-1;
    tm[end] = ne;
    nm = [nk; tm];
    #nm = [nk[1:end]; [nk[2:end]-1 ne]];
    if (nm[2,end]-nm[1,end]) < (na/2)
        nm[2,end-1] = nm[2,end];
        nm = nm[:,1:(end-1)];
    end
end
nm = Int.(nm);
nb = size(nm,2);

if nm[2,end]!=ne
    error("something wrong");
end
if maximum(nm[2,:]-nm[1,:])>ns
    error("something wrong");
end
if minimum(nm[2,:]-nm[1,:])<(ns/2-1)
    error("something wrong");
end

return nm,nb

end
