function facepartitionhdg(dmd,t,f,bcm,dim,elemtype,porder,nproc)

for i = 1:nproc
    println("face partition ", i)

    fi = f[:, dmd[i].elempart[:]]  # Ensure correct indexing in Julia
    dmd[i].bf = 0*fi;  # Initialize boundary flags    
    for j = 1:length(bcm)
        ind = findall(x -> x == j, fi)  # Find indices where fi equals j
        dmd[i].bf[ind] .= bcm[j]  # Assign bcm(j) to the corresponding indices
    end
        
    f2t,~ = mkf2e(t[:,dmd[i].elempart[:]],elemtype,dim);

    # only on [intelem bndelem]
    nelem = dmd[i].elempartpts;
    if length(nelem)>=2
        ne3 = sum(nelem[1:2]);
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
    dmd[i].f2t = f2t;    

    # fix f2t to ensure DOF consistency across interior faces
    # the global ID of the FIRST face must be smaller than that of the SECOND face        
    N = size(dmd[i].f2t, 2)
    for face in 1:N  # loop over each face
        e1 = dmd[i].f2t[1, face]
        e2 = dmd[i].f2t[3, face]
        if e2 > 0  # if the face is an interior face
            g1 = dmd[i].elempart[e1]
            g2 = dmd[i].elempart[e2]
            if g2 < g1
                # ensure g2 > g1
                tm = dmd[i].f2t[3:4, face]
                dmd[i].f2t[3:4, face] = dmd[i].f2t[1:2, face]
                dmd[i].f2t[1:2, face] = tm
            end
            e1 = dmd[i].f2t[1, face]
            e2 = dmd[i].f2t[3, face]
            g1 = dmd[i].elempart[e1]
            g2 = dmd[i].elempart[e2]
            if g2 <= g1
                error("something wrong")
            end
        end
    end

    dmd[i].facecon, dmd[i].elemcon = faceconnectivity2(t[:,dmd[i].elempart[:]],dmd[i].f2t,dim,elemtype,porder);
end

return dmd

end
