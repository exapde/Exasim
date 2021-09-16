function  exasim(pde,mesh)

if isa(pde, Array)
    nmodels = length(pde);      
else
    nmodels = 1;   
end

if nmodels==1
    # search compilers and set options
    pde = Gencode.setcompilers(pde);

    # generate input files and store them in datain folder
    pde, mesh, master, dmd = Preprocessing.preprocessing(pde,mesh);

    # generate source codes and store them in app folder
    Gencode.gencode(pde);

    if pde.usecmake==1    
        runstr = Main.cmakecompile(pde, nmodels); # use cmake to compile source codes 
        compilerstr = "";    
    else
        # compile source codes to build an executable file and store it in app folder
        compilerstr = Gencode.compilecode(pde);

        # run executable file to compute solution and store it in dataout folder
        runstr = Gencode.runcode(pde, nmodels);
    end

    # get solution from output files in dataout folder
    sol = Postprocessing.fetchsolution(pde,master,dmd,"dataout");
else
    master = Array{Any, 1}(undef, nmodels);
    dmd = Array{Any, 1}(undef, nmodels);
    sol = Array{Any, 1}(undef, nmodels);

    # preprocess and generate code for all PDE models
    for m = 1:nmodels    
        pde[m],mesh[m],master[m],dmd[m] = generatecode(pde[m],mesh[m]);
    end        
    gencodeall(nmodels);

    if pde.usecmake==1    
        runstr = Main.cmakecompile(pde[1], nmodels); # use cmake to compile source codes 
        compilerstr = "";    
    else
        # # compile source codes to build an executable file and store it in app folder
        compilerstr = compilecode(pde[1]);

        # run executable file to compute solution and store it in dataout folder
        runstr = runcode(pde[1],nmodels);
    end

    # get solution from output files in dataout folder
    for m = 1:nmodels        
        sol[m] = fetchsolution(pde[m], master[m], dmd[m], "dataout" * string(m));
    end    
end

return sol,pde,mesh,master,dmd,compilerstr,runstr

end
