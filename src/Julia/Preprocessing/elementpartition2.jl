function elementpartition2(dmd,t,t2t,nproc,metis)

nve,ne = size(t);

if nproc==1
    i = 1;
    dmd[i].elempart = reshape(1:ne,1,ne);
    dmd[i].elempartpts = reshape([ne],1,1);
    dmd[i].nbsd = reshape([0],1,1);
    #dmd[i].elem2cpu = reshape([0],1,1);
    dmd[i].elemrecv = reshape([0],1,1);
    dmd[i].elemsend = reshape([0],1,1);
    dmd[i].elemrecvpts = reshape([0],1,1);
    dmd[i].elemsendpts = reshape([0],1,1);
    return dmd;
end

elem2cpu,~ = partition(t,ne,nproc,metis);

for i = 1:nproc
    intelem = findall(elem2cpu[:] .== (i-1)); # elements in subdomain i
    elem = neighboringelements(t2t, intelem); # all elements connected to elements in subdomain i
    extelem = sort(setdiff(elem,intelem)); # exterior elements

    elem = neighboringelements(t2t, extelem); # all elements connected to exterior elements
    bndelem = intersect(elem,intelem);  # boundary elements in subdmain i
    outelem = sort(setdiff(elem,vcat(intelem, extelem))); #  elements outside subdmain i

    tmp = [setdiff(intelem,bndelem); bndelem; extelem; outelem];
    dmd[i].elempart = reshape(tmp,1,length(tmp)); # partitions of elements
    dmd[i].elempartpts = [length(intelem)-length(bndelem) length(bndelem) length(extelem) length(outelem)];
    #dmd[i].elem2cpu = elem2cpu[dmd[i].elempart];
    nelem = dmd[i].elempartpts;

    recvelem = [extelem[:]; outelem[:]]; # elements received from neighbors
    tmp = hcat(elem2cpu[recvelem].+1, Vector((sum(nelem[1:2])+1):sum(nelem)));
    dmd[i].elemrecv = hcat(tmp, recvelem);
    ind = sortperm(dmd[i].elemrecv[:,1]);
    dmd[i].elemrecv = dmd[i].elemrecv[ind,:];
    dmd[i].nbsd = (sort(unique(dmd[i].elemrecv[:,1])))'; # neighboring subdomains
end

# store elements sent to neighboring subdomains to assemble the linear system
for k = 1:nproc
    dmd[k].elemsend = zeros(Int,0,1);
end
for i = 1:nproc
    for j = 1:length(dmd[i].nbsd)
        # cpu k sends information to cpu i
        k = dmd[i].nbsd[j];
        ii = findall(dmd[i].elemrecv[:,1] .== k);
        tm = dmd[i].elemrecv[ii,:];
        tm[:,1] .= i;
        tm[:,2] = indexin(tm[:,3],dmd[k].elempart[:]);
        if length(dmd[k].elemsend)==0
            dmd[k].elemsend = tm;
        else
            dmd[k].elemsend = [dmd[k].elemsend; tm];
        end
    end
end

for i = 1:nproc
    dmd[i].elemsendpts = zeros(Int,1, length(dmd[i].nbsd));
    dmd[i].elemrecvpts = zeros(Int,1, length(dmd[i].nbsd));
    for j = 1:length(dmd[i].nbsd)
        dmd[i].elemsendpts[j] = length(findall(dmd[i].elemsend[:,1] .== dmd[i].nbsd[j]));
        dmd[i].elemrecvpts[j] = length(findall(dmd[i].elemrecv[:,1] .== dmd[i].nbsd[j]));
    end
    dmd[i].elemsend = reshape(dmd[i].elemsend[:,2],1,size(dmd[i].elemsend,1));
    dmd[i].elemrecv = reshape(dmd[i].elemrecv[:,2],1,size(dmd[i].elemrecv,1));
end

return dmd; # elem2cpu; #elempart,elempartall

end
