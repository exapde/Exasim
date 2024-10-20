function meshpartitionhdg(dmd, t,f,t2t,bcm,dim,elemtype,porder,nproc,metis)

display("run elementpartition...");  
dmd = elementpartitionhdg(dmd, t,t2t,nproc,metis);

display("run facepartition...");  
dmd = facepartitionhdg(dmd,t,f,bcm,dim,elemtype,porder,nproc);

return dmd

end


