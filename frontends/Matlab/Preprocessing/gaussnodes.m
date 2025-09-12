function [xgauss, wgauss] = gaussnodes(pgauss,dim,elemtype)

fileID = fopen("gaussnodes.bin",'r');
tmp = fread(fileID,'double');
fclose(fileID);

ndims = tmp(1);
k1 = 2; k2 = (k1-1)+(ndims);
narrays = tmp(k1:k2);
k1 = k2+1; k2 = (k1-1) + prod(narrays);
sz1 = reshape(tmp(k1:k2), narrays(:)');
k1 = k2+1; k2 = (k1-1) + prod(narrays);
sz2 = reshape(tmp(k1:k2), narrays(:)');
sz = sz1.*sz2;
lz = cumsum([0; sz(:)]);
lz1 = reshape(lz(1:end-1),narrays(:)');
lz2 = reshape(lz(2:end),narrays(:)');
      
e = elemtype+1;

i = 1;
m1 = k2+1 + lz1(i,e,pgauss,dim);
m2 = m1-1 + lz2(i,e,pgauss,dim)-lz1(i,e,pgauss,dim);
xgauss = reshape(tmp(m1:m2),sz1(i,e,pgauss,dim),sz2(i,e,pgauss,dim));

i = 2;
m1 = k2+1 + lz1(i,e,pgauss,dim);
m2 = m1-1 + lz2(i,e,pgauss,dim)-lz1(i,e,pgauss,dim);
wgauss = reshape(tmp(m1:m2),sz1(i,e,pgauss,dim),sz2(i,e,pgauss,dim));


