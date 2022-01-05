def fetchresidual(pde):

    if isinstance(pde, list):
        nmodels = len(pde);      
    else:
        nmodels = 1
    
    res = None;
    if nmodels == 1:        
        if pde['saveResNorm']:
            fn = "dataout/out_residualnorms0.bin";
            tm = fromfile(open(fn, "r"), dtype=float64);
            ne = int(round(size(tm)/(4)));            
            tm = reshape(tm,[4,ne],'F');                
            res = tm.transpose();

    else:        
        res = [None] * nmodels

        # get resiudal from output files in dataout folder
        for m in range(0, nmodels):
            if pde[m]['saveResNorm']:
                fn = "dataout/out_residualnorms" + str(m) + ".bin";
                tm = fromfile(open(fn, "r"), dtype=float64);
                ne = int(round(size(tm)/(4)));            
                tm = reshape(tm,[4,ne],'F');                
                res[m] = tm.transpose();

    return res

