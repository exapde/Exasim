function runstr = runcode(pde, numpde)

disp("Run C++ Exasim code ...")
cdir = pwd();
cd(char(pde.buildpath));

pdenum = " " + num2str(numpde) + " ";
DataPath = pde.buildpath;

mystr = pdenum;
if numpde==1
    mystr = mystr + DataPath + "/datain/ " + DataPath + "/dataout/out";
else    
    for i = 1:numpde
        mystr = mystr + DataPath + "/datain" + num2str(i) + "/ " + DataPath + "/dataout" + num2str(i) + "/out";
        mystr = mystr + " ";
    end    
end

mpirun = pde.mpirun;
if pde.platform == "cpu"        
    if pde.mpiprocs==1        
        exec = "./cpuEXASIM ";        
        runstr = "!" + exec + mystr;
    else        
        exec = " ./cpumpiEXASIM ";
        runstr = "!" + mpirun + " -np " + string(pde.mpiprocs) + exec + mystr;
    end        
elseif pde.platform == "gpu"
    if pde.mpiprocs==1        
        exec = "./gpuEXASIM ";     
        runstr = "!" + exec + mystr;
    else        
        exec = " ./gpumpiEXASIM ";
        runstr = "!" + mpirun + " -np " + string(pde.mpiprocs) + exec + mystr;
    end
end

tic
eval(char(runstr));
toc

cd(char(cdir));

end

