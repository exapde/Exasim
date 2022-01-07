function fetchresidual(pde)

if isa(pde, Array)
    nmodels = length(pde);      
else
    nmodels = 1;   
end

if nmodels==1
    if pde.saveResNorm
        fn = "dataout/out_residualnorms0.bin";
        res = reinterpret(Float64,read(fn));        
        ne = Int64(round(length(res)/4));
        res = reshape(res,(4,ne));                
        res = res';
    else
        res = [];
    end    
else
    res = Array{Any, 1}(undef, nmodels);
    # get residual from output files in dataout folder
    for m = 1:nmodels        
        if pde[m].saveResNorm
            fn = "dataout/out_residualnorms" * string(m-1) * ".bin";
            tm = reinterpret(Float64,read(fn));        
            ne = Int64(round(length(tm)/4));
            tm = reshape(tm,(4,ne));                
            res[m] = tm';
        end    
    end    
end

return res

end


