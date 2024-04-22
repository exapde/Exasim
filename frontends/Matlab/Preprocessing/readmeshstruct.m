function mesh = readmeshstruct(filemesh)
    
tm = readbin(filemesh);
sz = tm(1);
k1 = 2;
k2 = k1+(sz)-1;
mesh.nsize = tm(k1:k2);

k1 = k2+1;
k2 = k1+mesh.nsize(1)-1;
mesh.ndims = tm(k1:k2);

k1 = k2+1;
k2 = k1+mesh.nsize(2)-1;
mesh.facecon = tm(k1:k2);

k1 = k2+1;
k2 = k1+mesh.nsize(3)-1;
mesh.eblks = tm(k1:k2);

k1 = k2+1;
k2 = k1+mesh.nsize(4)-1;
mesh.fblks = tm(k1:k2);

k1 = k2+1;
k2 = k1+mesh.nsize(5)-1;
mesh.nbsd = tm(k1:k2);

k1 = k2+1;
k2 = k1+mesh.nsize(6)-1;
mesh.elemsend = tm(k1:k2);

k1 = k2+1;
k2 = k1+mesh.nsize(7)-1;
mesh.elemrecv = tm(k1:k2);

k1 = k2+1;
k2 = k1+mesh.nsize(8)-1;
mesh.elemsendpts = tm(k1:k2);

k1 = k2+1;
k2 = k1+mesh.nsize(9)-1;
mesh.elemrecvpts = tm(k1:k2);

k1 = k2+1;
k2 = k1+mesh.nsize(10)-1;
mesh.elempart = tm(k1:k2);

k1 = k2+1;
k2 = k1+mesh.nsize(11)-1;
mesh.elempartpts = tm(k1:k2);

k1 = k2+1;
k2 = k1+mesh.nsize(12)-1;
mesh.cgelcon = tm(k1:k2);

k1 = k2+1;
k2 = k1+mesh.nsize(13)-1;
mesh.rowent2elem = tm(k1:k2);

k1 = k2+1;
k2 = k1+mesh.nsize(14)-1;
mesh.cgent2dgent = tm(k1:k2);

k1 = k2+1;
k2 = k1+mesh.nsize(15)-1;
mesh.colent2elem = tm(k1:k2);

k1 = k2+1;
k2 = k1+mesh.nsize(16)-1;
mesh.rowe2f1 = tm(k1:k2);

k1 = k2+1;
k2 = k1+mesh.nsize(17)-1;
mesh.cole2f1 = tm(k1:k2);

k1 = k2+1;
k2 = k1+mesh.nsize(18)-1;
mesh.ent2ind1 = tm(k1:k2);

k1 = k2+1;
k2 = k1+mesh.nsize(19)-1;
mesh.rowe2f2 = tm(k1:k2);

k1 = k2+1;
k2 = k1+mesh.nsize(20)-1;
mesh.cole2f2 = tm(k1:k2);

k1 = k2+1;
k2 = k1+mesh.nsize(21)-1;
mesh.ent2ind2 = tm(k1:k2);


