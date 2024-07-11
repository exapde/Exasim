function masternodes(porder::Int,dim::Int,elemtype::Int)

d0 = pwd();
ii = findlast("Exasim", d0);

# filename = string(d0[1:ii[end]],"Preprocessing/masternodes.mat");
# vars = MAT.matread(filename);
#
# pelem1 = vars["pelem"][elemtype+1,porder,dim];
# telem1 = vars["telem"][elemtype+1,porder,dim];
# pface1 = vars["pface"][elemtype+1,porder,dim];
# tface1 = vars["tface"][elemtype+1,porder,dim];
# perm1 = vars["perm"][elemtype+1,porder,dim];
#
# telem1 = convert(Array{Int}, telem1);
# tface1 = convert(Array{Int}, tface1);
# perm1 = convert(Array{Int}, perm1);

fn = string(d0[1:ii[end]],"/frontends/Julia/Preprocessing/masternodes.bin");
tmp = reinterpret(Float64,read(fn));

ndims = Int(tmp[1]);
k1 = 2; k2 = (k1-1)+(ndims);
narrays = Int.(tmp[k1:k2]);
k1 = k2+1; k2 = (k1-1) + prod(narrays);
sz1 = reshape(Int.(tmp[k1:k2]), (narrays...,));
k1 = k2+1; k2 = (k1-1) + prod(narrays);
sz2 = reshape(Int.(tmp[k1:k2]), (narrays...,));
sz = sz1.*sz2;
lz = cumsum([0; sz[:]]);
lz1 = reshape(lz[1:end-1],(narrays...,));
lz2 = reshape(lz[2:end],(narrays...,));

e = elemtype+1;
i = 1;
m1 = k2+1 + lz1[i,e,porder,dim];
m2 = m1-1 + lz2[i,e,porder,dim]-lz1[i,e,porder,dim];
pelem = reshape(tmp[m1:m2],sz1[i,e,porder,dim],sz2[i,e,porder,dim]);

i = 2;
m1 = k2+1 + lz1[i,e,porder,dim];
m2 = m1-1 + lz2[i,e,porder,dim]-lz1[i,e,porder,dim];
telem = reshape(tmp[m1:m2],sz1[i,e,porder,dim],sz2[i,e,porder,dim]);

i = 3;
m1 = k2+1 + lz1[i,e,porder,dim];
m2 = m1-1 + lz2[i,e,porder,dim]-lz1[i,e,porder,dim];
pface = reshape(tmp[m1:m2],sz1[i,e,porder,dim],sz2[i,e,porder,dim]);

i = 4;
m1 = k2+1 + lz1[i,e,porder,dim];
m2 = m1-1 + lz2[i,e,porder,dim]-lz1[i,e,porder,dim];
tface = reshape(tmp[m1:m2],sz1[i,e,porder,dim],sz2[i,e,porder,dim]);

i = 5;
m1 = k2+1 + lz1[i,e,porder,dim];
m2 = m1-1 + lz2[i,e,porder,dim]-lz1[i,e,porder,dim];
perm = reshape(tmp[m1:m2],sz1[i,e,porder,dim],sz2[i,e,porder,dim]);

telem = convert(Array{Int}, telem);
tface = convert(Array{Int}, tface);
perm = convert(Array{Int}, perm);

# display(pelem - pelem1)
# display(telem - telem1)
# display(pface - pface1)
# display(tface - tface1)
# display(perm - perm1)

return pelem, telem, pface, tface, perm;

end
