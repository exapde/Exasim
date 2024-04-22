function dmd = meshpartitionhdg(t,f,t2t,bcm,dim,elemtype,porder,nproc,metis)

disp('run elementpartition...');  
dmd = elementpartitionhdg(t,t2t,nproc,metis);

disp('run facepartition...');  
dmd = facepartitionhdg(dmd,t,f,bcm,dim,elemtype,porder,nproc);



