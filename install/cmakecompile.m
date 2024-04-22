function cmakecompile(pde, numpde)

cdir = pwd();
DataPath = pde.buildpath;
cd(char(pde.buildpath));

pdenum = " " + num2str(numpde) + " ";

if exist("cpuEXASIM", "file")
    delete(char("cpuEXASIM"));
end
if exist("cpumpiEXASIM", "file")
    delete(char("cpumpiEXASIM"));
end
if exist("gpuEXASIM", "file")
    delete(char("gpuEXASIM"));
end
if exist("gpumpiEXASIM", "file")
    delete(char("gpumpiEXASIM"));
end
if exist("cpuApp.a", "file")
    delete(char("cpuApp.a"));
end
if exist("gpuApp.a", "file")
    delete(char("gpuApp.a"));
end

if pde.platform == "gpu"
  eval("!cmake -D EXASIM_APP=ON -D EXASIM_CORES=ON -D EXASIM_EXECUTABLES=ON -D EXASIM_CUDA=ON ../install");
else
  eval("!cmake -D EXASIM_APP=ON -D EXASIM_CORES=ON -D EXASIM_EXECUTABLES=ON ../install");
end
eval("!cmake --build .");

tic
mpirun = pde.mpirun;
if pde.platform == "cpu"        
    if pde.mpiprocs==1        
        str = "!./cpuEXASIM " + pdenum + DataPath + "/datain/ " + DataPath + "/dataout/out";        
    else        
        str = "!" + mpirun + " -np " + string(pde.mpiprocs) + " ./cpumpiEXASIM " + pdenum + DataPath + "/datain/ " + DataPath + "/dataout/out";       
    end    
    eval(char(str));
elseif pde.platform == "gpu"
    if pde.mpiprocs==1        
        str = "!./gpuEXASIM " + pdenum + DataPath + "/datain/ " + DataPath + "/dataout/out";        
    else        
        str = "!" + mpirun + " -np " + string(pde.mpiprocs) + " ./app/gpumpiEXASIM " + pdenum + DataPath + "/datain/ " + DataPath + "/dataout/out";       
    end
    eval(char(str));
end

cd(char(cdir));

toc