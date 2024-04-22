import Preprocessing, Gencode
from fetchsolution import fetchsolution
from generatecode import generatecode

def exasim(pde,mesh):

    if isinstance(pde, list):
        nmodels = len(pde);      
    else:
        nmodels = 1
    
    res = None;
    if nmodels == 1:
        # search compilers and set options
        pde = Gencode.setcompilers(pde);

        # generate input files and store them in datain folder
        pde, mesh, master, dmd = Preprocessing.preprocessing(pde,mesh);

        # generate source codes and store them in app folder
        Gencode.gencode(pde);

        # compile source codes to build an executable file and store it in app folder
        compilerstr = Gencode.compilecode(pde);

        # run executable file to compute solution and store it in dataout folder
        runstr = Gencode.runcode(pde, nmodels);

        # get solution from output files in dataout folder
        pde['vistime'] = [];
        sol = fetchsolution(pde,master,dmd,"dataout");
        
        if pde['saveResNorm']:
            fn = "dataout/out_residualnorms0.bin";
            tm = fromfile(open(fn, "r"), dtype=float64);
            ne = int(round(size(tm)/(4)));            
            tm = reshape(tm,[4,ne],'F');                
            res = tm.transpose();

    else:        
        master = [None] * nmodels
        dmd = [None] * nmodels
        sol = [None] * nmodels
        res = [None] * nmodels

        # preprocess and generate code for all PDE models
        for m in range(0, nmodels):    
            pde[m],mesh[m],master[m],dmd[m] = generatecode(pde[m],mesh[m])[0:4];
        
        Gencode.gencodeall(nmodels);

        # compile source codes to build an executable file and store it in app folder
        compilerstr = Gencode.compilecode(pde[0]);

        # run executable file to compute solution and store it in dataout folder
        runstr = Gencode.runcode(pde[0],nmodels);

        # get solution from output files in dataout folder
        for m in range(0, nmodels):
            sol[m] = fetchsolution(pde[m],master[m],dmd[m], "dataout" + str(m+1));        

            if pde[m]['saveResNorm']:
                fn = "dataout/out_residualnorms" + str(m) + ".bin";
                tm = fromfile(open(fn, "r"), dtype=float64);
                ne = int(round(size(tm)/(4)));            
                tm = reshape(tm,[4,ne],'F');                
                res[m] = tm.transpose();

    return sol,pde,mesh,master,dmd,compilerstr,runstr,res
