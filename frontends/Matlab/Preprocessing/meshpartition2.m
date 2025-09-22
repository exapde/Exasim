function dmd = meshpartition2(t,f,t2t,bcm,dim,elemtype,porder,nproc,metis,Cxxpreprocessing)

disp('run elementpartition...');  
dmd = elementpartition2(t,t2t,nproc,metis);

if Cxxpreprocessing==0
    disp('run facepartition...');  
    dmd = facepartition2(dmd,t,f,bcm,dim,elemtype,porder,nproc);
else
    for i = 1:nproc   
        fi = f(:,dmd{i}.elempart); 
        dmd{i}.bf = 0*f(:,dmd{i}.elempart);
        for j=1:length(bcm)
          ind = fi==j;
          dmd{i}.bf(ind) = bcm(j);
        end      
    end    
end

