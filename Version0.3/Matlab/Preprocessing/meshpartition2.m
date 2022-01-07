function dmd = meshpartition2(t,f,t2t,bcm,dim,elemtype,porder,nproc,metis)

disp('run elementpartition...');  
dmd = elementpartition2(t,t2t,nproc,metis);

disp('run facepartition...');  
dmd = facepartition2(dmd,t,f,bcm,dim,elemtype,porder,nproc);




