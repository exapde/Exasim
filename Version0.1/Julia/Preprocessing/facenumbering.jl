function  facenumbering(p,t,elemtype,bndexpr,prdexpr)

dim = size(p,1);
ne = size(t,2);
face = getelemface(dim,elemtype);
nvf,nfe = size(face);

t2fl = reshape(t[face[:],:], nvf, nfe*ne);
pf = reshape(p[:,t2fl[:]], dim, nvf, nfe*ne);
pf = reshape(sum(pf,dims=2)/nvf,dim,nfe,ne);

# interior faces are zero
f = zeros(Int,nfe,ne);

# face-to-element connectivity
#f2t = mkf2t(t,elemtype,dim);
f2t,t2t = mkf2e(t,elemtype,dim);

# find elements on the domain boundary
ind = findall(f2t[3,:].==0);
for i = 1:length(ind) # for each element on the domain boundary
    e = f2t[1,ind[i]]; # element e
    l = f2t[2,ind[i]]; # local face index
    for k = 1:length(bndexpr) # for each boundary expression
        a = bndexpr[k](pf[:,l,e]); # evaluate the boundary expression
        if a[1] # check if element e belong to this boundary
            f[l,e] = k; # then set f(l,e) to k
            break;
        end
    end
end
# # regular boundary faces are positive
# for i = 1:length(bndexpr)
#     tm = bndexpr[i](pf);
#     f[tm .== 1] .= i;
# end

nprd = size(prdexpr,1);
if nprd>0
    f = f[:];
    pf = reshape(pf, dim, nfe*ne);
    tprd = copy(t);
    # periodic boundary faces are negative
    for i = 1:nprd
        i1 = findall(f .== prdexpr[i,1]);
        f[i1] = -f[i1];
        p1 = prdexpr[i,2](pf[:,i1[:]]);
        e1 = Int.(ceil.(Float64.(i1)./Float64(nfe)));      # 1st elements
        l1 = i1 .- (e1.-1).*nfe;   # 1st local faces

        i2 = findall(f .== prdexpr[i,3]);
        f[i2] = -f[i2];
        p2 = prdexpr[i,4](pf[:,i2[:]]);
        e2 = Int.(ceil.(Float64.(i2)./Float64(nfe)));      # 1st elements
        l2 = i2 .- (e2.-1).*nfe;   # 1st local faces

        # update t2t to connect periodic elements
        for j=1:length(e1)
            t2t[l2[j],e2[j]] = e1[j];
            t2t[l1[j],e1[j]] = e2[j];
        end

        p1 = reshape(p1, Int64(length(p1)/length(i1)), length(i1));
        p2 = reshape(p2, Int64(length(p2)/length(i2)), length(i2));
        in = xiny(p1,p2,0);

        v1 = t2fl[:,i1];
        v1 = sort(unique(v1[:]));
        p1 = prdexpr[i,2](p[:,v1]);

        v2 = t2fl[:,i2[in]];
        v2 = sort(unique(v2[:]));
        p2 = prdexpr[i,4](p[:,v2]);

        p1 = reshape(p1, Int64(length(p1)/length(v1)), length(v1));
        p2 = reshape(p2, Int64(length(p2)/length(v2)), length(v2));
        in = xiny(p1,p2,0);
        v2 = v2[in];

        for j=1:length(v1)
            tprd[tprd .== v2[j]] .= v1[j];
        end
    end
else
    #tprd = zeros(1,1);
    tprd = t;
end

f = reshape(f, nfe, ne);

return f, tprd, t2t

end
