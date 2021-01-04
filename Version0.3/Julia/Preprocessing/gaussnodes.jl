function gaussnodes(pgauss::Int,dim::Int,elemtype::Int)

d0 = pwd();
ii = findlast("Exasim", d0);

# filename = string(d0[1:ii[end]],"Preprocessing/gaussnodes.mat");
# vars = MAT.matread(filename);
#
# xgauss1 = vars["xgauss"][elemtype+1,pgauss,dim];
# wgauss1 = vars["wgauss"][elemtype+1,pgauss,dim];

#fn = string(d0[1:ii[end]],"Preprocessing/gaussnodes.bin");
fn = string(d0[1:ii[end]],"/Version0.1/Julia/Preprocessing/gaussnodes.bin");
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
m1 = k2+1 + lz1[i,e,pgauss,dim];
m2 = m1-1 + lz2[i,e,pgauss,dim]-lz1[i,e,pgauss,dim];
xgauss = reshape(tmp[m1:m2],sz1[i,e,pgauss,dim],sz2[i,e,pgauss,dim]);
i = 2;
m1 = k2+1 + lz1[i,e,pgauss,dim];
m2 = m1-1 + lz2[i,e,pgauss,dim]-lz1[i,e,pgauss,dim];
wgauss = reshape(tmp[m1:m2],sz1[i,e,pgauss,dim],sz2[i,e,pgauss,dim]);

# display(xgauss - xgauss1)
# display(wgauss - wgauss1)

return xgauss, wgauss;

end
