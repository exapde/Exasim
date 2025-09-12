function [pelem, telem, pface, tface, perm] = masternodes(porder,dim,elemtype)

fileID = fopen("masternodes.bin",'r');
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
m1 = k2+1 + lz1(i,e,porder,dim);
m2 = m1-1 + lz2(i,e,porder,dim)-lz1(i,e,porder,dim);
pelem = reshape(tmp(m1:m2),sz1(i,e,porder,dim),sz2(i,e,porder,dim));

i = 2;
m1 = k2+1 + lz1(i,e,porder,dim);
m2 = m1-1 + lz2(i,e,porder,dim)-lz1(i,e,porder,dim);
telem = reshape(tmp(m1:m2),sz1(i,e,porder,dim),sz2(i,e,porder,dim));

i = 3;
m1 = k2+1 + lz1(i,e,porder,dim);
m2 = m1-1 + lz2(i,e,porder,dim)-lz1(i,e,porder,dim);
pface = reshape(tmp(m1:m2),sz1(i,e,porder,dim),sz2(i,e,porder,dim));

i = 4;
m1 = k2+1 + lz1(i,e,porder,dim);
m2 = m1-1 + lz2(i,e,porder,dim)-lz1(i,e,porder,dim);
tface = reshape(tmp(m1:m2),sz1(i,e,porder,dim),sz2(i,e,porder,dim));

i = 5;
m1 = k2+1 + lz1(i,e,porder,dim);
m2 = m1-1 + lz2(i,e,porder,dim)-lz1(i,e,porder,dim);
perm = reshape(tmp(m1:m2),sz1(i,e,porder,dim),sz2(i,e,porder,dim));

perm = int64(perm);
tface = int64(tface);
telem = int64(telem);



