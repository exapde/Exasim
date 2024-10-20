function runcode(pde, numpde=1)

display("Run C++ Exasim code ...")

cdir = pwd(); 
cd(pde.buildpath);

buildpath = pde.buildpath;
DataPath = buildpath;
mpirun = pde.mpirun;
pdenum = string(numpde) * " ";

if pde.platform == "cpu"        
    if pde.mpiprocs==1        
        str = "./cpuEXASIM " * pdenum * DataPath * "/datain/ " * DataPath * "/dataout/out";        
    else        
        str = mpirun * " -np " * string(pde.mpiprocs) * " ./cpumpiEXASIM " * pdenum * DataPath * "/datain/ " * DataPath * "/dataout/out";       
    end        
elseif pde.platform == "gpu"
    if pde.mpiprocs==1        
        str = "./gpuEXASIM " * pdenum * DataPath * "/datain/ " * DataPath * "/dataout/out";        
    else        
        str = mpirun * " -np " * string(pde.mpiprocs) * " ./gpumpiEXASIM " * pdenum * DataPath * "/datain/ " * DataPath * "/dataout/out";       
    end    
end

run(string2cmd(str));

cd(cdir);

return str;

end
