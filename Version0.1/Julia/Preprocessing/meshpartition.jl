function meshpartition(dmd,p,t,f,t2t,tprd,elemtype,bcm,bndexpr,prdexpr,porder,nproc,metis)

# display("run facenumbering...");
# f, tprd = facenumbering(p,t,elemtype,bndexpr,prdexpr);

print("run elementpartition...\n");
if size(tprd,2)==size(t,1)
    dmd = elementpartition2(dmd,tprd,t2t,nproc,metis);
else
    dmd = elementpartition2(dmd,t,t2t,nproc,metis);
end

print("run facepartition...\n");
dmd = facepartition(dmd,p,t,tprd,f,bcm,elemtype,prdexpr,porder,nproc);

return dmd;

end
