function str = runcode(app)

disp("run code...");
mpirun = app.mpirun;
if app.platform == "cpu"    
    if app.mpiprocs==1
        str = "!./app/serial" + app.appname + " datain/ dataout/out";        
    else
        str = "!" + mpirun + " -np " + string(app.mpiprocs) + " ./app/mpi" + app.appname + " datain/ dataout/out";        
    end
    eval(char(str));
elseif app.platform == "gpu"
    if app.mpiprocs==1
        str = "!./app/gpu" + app.appname + " datain/ dataout/out";   
    else
        str = "!" + mpirun + " -np " + string(app.mpiprocs) + " ./app/gpumpi" + app.appname + " datain/ dataout/out";        
    end
    eval(char(str));
end

end

