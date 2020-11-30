function  mkv2t(t,ne)

sz = size(t);
if sz[1]==ne
    t = t';
end

ndof = maximum(t[:]);
re = zeros(Int,ndof,1); # store number of neighboring elements for each entity
for i = 1:ne
    k = t[:,i];
    k = k[k.>0];
    re[k] = re[k] .+ 1;
end
re = [0; cumsum(re,dims=1)];

ce = zeros(Int,re[end],1);  # store neighboring-element indices for each entity
in = ones(Int,ndof,1);
for i = 1:ne
    k = t[:,i];   # entities on element i
    k = k[k.>0];
    # re(k): pointer to the list of entities k
    ce[re[k]+in[k]] .= i;
    in[k] = in[k] .+ 1;
end

return re,ce

end
