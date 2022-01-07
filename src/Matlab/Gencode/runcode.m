function runstr = runcode(app, numpde)

if nargin<2
    numpde=1;
end

disp("run code...");
if numpde==1
    mystr = app.appname + " 1 datain/ dataout/out";
else
    mystr = app.appname + " " + num2str(numpde);
    for i = 1:numpde
        mystr = mystr + " datain" + num2str(i) + "/ " + "dataout" + num2str(i) + "/out";
    end    
end

mpirun = app.mpirun;
if app.platform == "cpu"    
    if app.mpiprocs==1
        runstr = "!./app/serial" + mystr;        
    else
        runstr = "!" + mpirun + " -np " + string(app.mpiprocs) + " ./app/mpi" + mystr;        
    end
    eval(char(runstr));
elseif app.platform == "gpu"
    if app.mpiprocs==1
        runstr = "!./app/gpu" + mystr;   
    else
        runstr = "!" + mpirun + " -np " + string(app.mpiprocs) + " ./app/gpumpi" + mystr;        
    end
    eval(char(runstr));
end

end

