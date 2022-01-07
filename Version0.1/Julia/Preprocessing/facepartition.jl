function facepartition(dmd,p,t,tprd,f,bcm,elemtype,prdexpr,porder,nproc)

nve,ne = size(t);
dim = size(p,1);
# if nve==2
#     dim=1;
# else
#     if elemtype==0 # tri/tet elements
#         dim=nve-1;
#     elseif elemtype==1 # quad/hex elements
#         dim=log2(nve);
#     end
# end

for i = 1:nproc
    fi = f[:,dmd[i].elempart[:]];
    if size(tprd,2)==ne
        #f2t = mkf2t(tprd[:,dmd[i].elempart[:]],elemtype,dim);
        f2t,~ = mkf2e(tprd[:,dmd[i].elempart[:]],elemtype,dim);
    else
        #f2t = mkf2t(t[:,dmd[i].elempart[:]],elemtype,dim);
        f2t,~ = mkf2e(t[:,dmd[i].elempart[:]],elemtype,dim);
    end

    # only on [intelem bndelem extelem]
    nelem = dmd[i].elempartpts;
    if length(nelem)>=3
        ne3 = sum(nelem[1:3]);
        ind = findall(f2t[1,:].<=ne3);
        f2t = f2t[:,ind];
    end

    # reorder so that boundary faces are last
    ina = findall(f2t[3,:].>0); # interior faces
    inb = findall(f2t[3,:].==0); # boundary faces
    #inc = sub2ind(size(fi), f2t[2,inb], f2t[1,inb]);
    inc = (f2t[1,inb].-1)*(size(fi,1)) .+ f2t[2,inb];
    fb = fi[inc]; # list of boundary indices

    fa = sort(unique(fb[:])); # boundary indices
    bcn = sort(unique(bcm[fa])); # a list of boundary conditions
    nbc = length(bcn);

    dmd[i].facepartpts = reshape([length(ina)],1,1);
    dmd[i].facepartbnd = reshape([0],1,1);
    ind = zeros(Int,1,length(fb));
    m = 1;
    for j=1:nbc # for each boundary condition bcn(j)
        bj = findall(bcm[:].==bcn[j]); # find all boundaries that have condition bcn(j)
        n = 0;
        for k = 1:length(bj) # for each boundary that has condition bcn(j)
            ii = findall(fb[:] .== bj[k]); # indices of the boundary bj(k)
            l = length(ii);
            n = n + l;
            ind[m:(m+l-1)] = ii;
            m = m + l;
        end
        dmd[i].facepartpts = [dmd[i].facepartpts n];
        dmd[i].facepartbnd = [dmd[i].facepartbnd bcn[j]];
    end
    # [interior faces, boundary faces]
    f2t = f2t[:,[ina; inb[ind[:]]]];

    dmd[i].facecon = faceconnectivity(p,t[:,dmd[i].elempart[:]],fi,f2t,elemtype,prdexpr,porder);
end

return dmd

end
