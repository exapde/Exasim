function  exasim(pde,mesh)

cdir = pwd(); ii = findlast("Exasim", cdir);
ExasimPath = cdir[1:ii[end]];
buildpath = ExasimPath * "/build";

if isa(pde, Array)
    nmodels = length(pde);      
else
    nmodels = 1;   
end
res = []
if nmodels==1
    # search compilers and set options
    pde = Gencode.setcompilers(pde);

    # generate input files and store them in datain folder
    pde, mesh, master, dmd = Preprocessing.preprocessing(pde,mesh);

    # generate source codes and store them in app folder
    Gencode.gencode(pde);

    compilerstr = Main.cmakecompile(pde); # use cmake to compile source codes 
    #compilerstr = Main.compilepdemodel(pde); # use cmake to compile source codes 

    runstr = Gencode.runcode(pde, 1);

    # get solution from output files in dataout folder
    sol = Postprocessing.fetchsolution(pde,master,dmd, buildpath * "/dataout");
    if pde.saveResNorm == 1
        fn = "dataout/out_residualnorms0.bin";
        res = reinterpret(Float64,read(fn));        
        ne = Int64(round(length(res)/4));
        res = reshape(res,(4,ne));                
        res = res';
    end    
else
    master = Array{Any, 1}(undef, nmodels);
    dmd = Array{Any, 1}(undef, nmodels);
    sol = Array{Any, 1}(undef, nmodels);
    res = Array{Any, 1}(undef, nmodels);

    # preprocess and generate code for all PDE models
    for m = 1:nmodels    
        pde[m],mesh[m],master[m],dmd[m] = generatecode(pde[m],mesh[m]);
    end        
    gencodeall(nmodels);

    if pde.usecmake==1    
        compilerstr = Main.cmakecompile(pde[1]); # use cmake to compile source codes 
        #compilerstr = "";    
        #compilerstr = Main.compilepdemodel(pde[1]); # use cmake to compile source codes 
        runstr = runcode(pde[1],nmodels);
    else
        # # compile source codes to build an executable file and store it in app folder
        compilerstr = compilecode(pde[1]);

        # run executable file to compute solution and store it in dataout folder
        runstr = runcode(pde[1],nmodels);
    end

    # get solution from output files in dataout folder
    for m = 1:nmodels        
        sol[m] = Postprocessing.fetchsolution(pde[m], master[m], dmd[m], "dataout" * string(m));
        #sol[m],~,~ = Postprocessing.getsolutions(pde[m], dmd[m]);

        if pde[m].saveResNorm
            fn = "dataout/out_residualnorms" * string(m-1) * ".bin";
            tm = reinterpret(Float64,read(fn));        
            ne = Int64(round(length(tm)/4));
            tm = reshape(tm,(4,ne));                
            res[m] = tm';
        end    
    end    
end

return sol,pde,mesh,master,dmd,compilerstr,runstr,res

end

