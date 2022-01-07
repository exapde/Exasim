function runcode(app, numpde)

display("run code...");

if numpde==1
    mystr = app.appname * " 1 datain/ dataout/out";
else
    mystr = app.appname * " " * string(numpde);
    for i = 1:numpde
        mystr = mystr * " datain" * string(i) * "/ " * "dataout" * string(i) * "/out";
    end    
end

mpirun = app.mpirun;
if app.platform == "cpu"
    if app.mpiprocs==1
        runstr = "./app/serial" * mystr; 
    else
        runstr = mpirun * " -np " * string(app.mpiprocs) * " ./app/mpi" * mystr;        
    end
    run(string2cmd(runstr));
elseif app.platform == "gpu"
    if app.mpiprocs==1
        runstr = "./app/gpu" * mystr;   
    else
        runstr = mpirun * " -np " * string(app.mpiprocs) * " ./app/gpumpi" * mystr;                
    end
    run(string2cmd(runstr));
else
    error("app.platform must be either cpu or gpu");
end

return runstr

end
