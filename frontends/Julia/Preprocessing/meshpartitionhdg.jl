function meshpartitionhdg(dmd, t,f,t2t,bcm,dim,elemtype,porder,nproc,metis,Cxxpreprocessing)

display("run elementpartition...");  
dmd = elementpartitionhdg(dmd, t,t2t,nproc,metis);

if Cxxpreprocessing == 0    
    display("run facepartition...");  
    dmd = facepartitionhdg(dmd,t,f,bcm,dim,elemtype,porder,nproc);
else
    for i = 1:nproc   
        fi = f[:, dmd[i].elempart[:]]  # Ensure correct indexing in Julia
        dmd[i].bf = 0*fi;  # Initialize boundary flags    
        for j = 1:length(bcm)
            ind = findall(x -> x == j, fi)  # Find indices where fi equals j
            dmd[i].bf[ind] .= bcm[j]  # Assign bcm(j) to the corresponding indices
        end
    end
end

return dmd

end


