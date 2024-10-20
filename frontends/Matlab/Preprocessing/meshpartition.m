function dmd = meshpartition(p,t,f,t2t,tprd,elemtype,bcm,bndexpr,prdexpr,porder,nproc,metis)

% disp('run facenumbering...');  
% [f, tprd] = facenumbering(p,t,elemtype,bndexpr,prdexpr);

disp('run elementpartition...');  
if size(tprd,2)==size(t,1)
    dmd = elementpartition2(tprd,t2t,nproc,metis);
else
    dmd = elementpartition2(t,t2t,nproc,metis);
end

disp('run facepartition...');  
dmd = facepartition(dmd,p,t,tprd,f,bcm,elemtype,prdexpr,porder,nproc);



