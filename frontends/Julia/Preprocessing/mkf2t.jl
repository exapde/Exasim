function mkf2t(t,elemtype,dim)

ne = size(t,2);
face = getelemface(dim,elemtype);
nvf,nfe = size(face);

N = nfe*ne;
sfaces = reshape(t[face[:],:],nvf, N);
sfaces = sort(sfaces, dims=1);
f = 0*sfaces;
f[:,1] = sfaces[:,1];
f2t = zeros(Int,4,N);
f2t[1,1] = 1;
f2t[2,1] = 1;

k = 1;
for i = 2:N
    l = rem(i,nfe); # local face l
    if (l == 0)
        l = nfe;
    end
    e = (i-l)/nfe + 1; # element e

    diff = sum(abs.(sfaces[:,i] .- f[:,1:k]),dims=1);
    dmin = minimum(diff);
    if dmin==0 # match
        j = findall(diff[:] .== dmin);
        f2t[3,j] .= e;
        f2t[4,j] .= l;
    else # not mach
        k = k + 1;
        # add local face to the list
        f[:,k] = sfaces[:,i];
        f2t[1,k] = e;
        f2t[2,k] = l;
    end
end
f2t = f2t[:,1:k];

return f2t

end
