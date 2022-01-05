function meshpartition2(dmd,t,f,t2t,bcm,dim,elemtype,porder,nproc,metis)

print("run elementpartition...\n");
dmd = elementpartition2(dmd,t,t2t,nproc,metis);

print("run facepartition...\n");
dmd = facepartition2(dmd,t,f,bcm,dim,elemtype,porder,nproc);

return dmd;

end
