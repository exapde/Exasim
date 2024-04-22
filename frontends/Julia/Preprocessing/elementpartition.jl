function elementpartition(dmd,t,elemtype,nproc,metis)

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

ne = size(t,2);
elem2cpu,~ = partition(t,ne,nproc,metis);

if nve==2
    dim=1;
else
    if elemtype==0 # tri/tet elements
        dim=nve-1;
    elseif elemtype==1 # quad/hex elements
        dim=log2(nve);
    end
end
face = getelemface(dim,elemtype);
nvf,nfe = size(face);

display("run mkv2t...");
re,ce = mkv2t(t,ne);

for i = 1:nproc
    intelem = findall(elem2cpu[:] .== (i-1)); # elements in subdomain i
    elem = node2elem(t[:,intelem],re,ce);
    extelem = sort(setdiff(elem,intelem)); # exterior elements

    # fix extelem
    n1 = length(intelem);
    n2 = length(extelem);
    t1 = reshape(t[face,intelem],nvf, nfe*n1);
    t2 = reshape(t[face,extelem],nvf, nfe, n2);
    t1 = sort(t1,dims=1);
    t2 = sort(t2,dims=1);
    match = zeros(Int,1,n2);
    for j=1:n2
        for k=1:nfe
            if minimum(sum(abs.(reshape(t2[:,k,j],nvf,1) .- t1),dims=1))==0
                match[j] = 1;
                break;
            end
        end
    end
    ind = findall(match[:].==1);
    extelem = extelem[ind];

    elem = node2elem(t[:,extelem],re,ce); # all elements connected to exterior elements
    bndelem = intersect(elem,intelem);  # boundary elements in subdmain i
    outelem = sort(setdiff(elem,vcat(intelem, extelem))); #  elements outside subdmain i

    # fix outelem
    n1 = length([intelem; extelem]);
    n2 = length(outelem);
    t1 = reshape(t[face[:],[intelem; extelem]], nvf, nfe*n1);
    t2 = reshape(t[face[:],outelem], nvf, nfe, n2);
    t1 = sort(t1,dims=1);
    t2 = sort(t2,dims=1);
    match = zeros(Int,1,n2);
    for j=1:n2
        for k=1:nfe
            if minimum(sum(abs.(t2[:,k,j] .- t1),dims=1))==0
                match[j] = 1;
                break;
            end
        end
    end
    ind = findall(match[:].==1);
    outelem = outelem[ind];

    tmp = [setdiff(intelem,bndelem); bndelem; extelem; outelem];
    dmd[i].elempart = reshape(tmp,1,length(tmp)); # partitions of elements
    dmd[i].elempartpts = [length(intelem)-length(bndelem) length(bndelem) length(extelem) length(outelem)];
    #dmd[i].elem2cpu = elem2cpu[dmd[i].elempart];
    nelem = dmd[i].elempartpts;

    # recvelem = [extelem[:]; outelem[:]]; # elements received from neighbors
    # ind = xiny(recvelem[:], dmd{i}.elempart[:], 1);
    # elem2cpu = dmd[i].elem2cpu[ind[:]];
    # tmp = hcat(elem2cpu+1, Vector((sum(nelem[1:2])+1):sum(nelem)));
    # dmd[i].elemrecv = hcat(tmp, recvelem);
    # ind = sortperm(elem2cpu[:]);
    # dmd[i].elemrecv = dmd[i].elemrecv[ind,:];
    # dmd[i].nbsd = (sort(unique(dmd[i].elemrecv[:,1])))'; # neighboring subdomains

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
