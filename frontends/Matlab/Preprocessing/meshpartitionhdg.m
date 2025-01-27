function dmd = meshpartitionhdg(t,f,t2t,bcm,dim,elemtype,porder,coupledinterface,nproc,metis)

disp('run elementpartition...');  
dmd = elementpartitionhdg(t,t2t,f,coupledinterface,nproc,metis);

disp('run facepartition...');  
dmd = facepartitionhdg(dmd,t,f,bcm,dim,elemtype,porder,nproc);



