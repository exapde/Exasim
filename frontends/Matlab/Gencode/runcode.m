function runstr = runcode(pde, numpde)

disp("Run C++ Exasim code ...")
cdir = pwd();
cd(char(pde.buildpath));

pdenum = " " + num2str(numpde) + " ";
DataPath = pde.buildpath;

mpirun = pde.mpirun;
if pde.platform == "cpu"        
    if pde.mpiprocs==1        
        exec = "./cpuEXASIM ";        
        runstr = "!" + exec + pdenum + DataPath + "/datain/ " + DataPath + "/dataout/out";        
    else        
        exec = " ./cpumpiEXASIM ";
        runstr = "!" + mpirun + " -np " + string(pde.mpiprocs) + exec + pdenum + DataPath + "/datain/ " + DataPath + "/dataout/out";       
    end        
elseif pde.platform == "gpu"
    if pde.mpiprocs==1        
        exec = "./gpuEXASIM ";     
        runstr = "!" + exec + pdenum + DataPath + "/datain/ " + DataPath + "/dataout/out";        
    else        
        exec = " ./gpumpiEXASIM ";
        runstr = "!" + mpirun + " -np " + string(pde.mpiprocs) + exec + pdenum + DataPath + "/datain/ " + DataPath + "/dataout/out";       
    end
end

tic
eval(char(runstr));
toc

cd(char(cdir));

end

