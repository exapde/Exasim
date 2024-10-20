function comstr = cmakecompile(pde)

disp("Compile C++ Exasim code using cmake...")

cdir = pwd();
cd(char(pde.buildpath));

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

if pde.mpiprocs==1 
  if pde.platform == "gpu"
    comstr = "!cmake -D EXASIM_NOMPI=ON -D EXASIM_MPI=OFF -D EXASIM_CUDA=ON ../install";    
  else
    comstr = "!cmake -D EXASIM_NOMPI=ON -D EXASIM_MPI=OFF -D EXASIM_CUDA=OFF ../install";
  end
else
  if pde.platform == "gpu"
    comstr = "!cmake -D EXASIM_NOMPI=OFF -D EXASIM_MPI=ON -D EXASIM_CUDA=ON ../install";  
  else
    comstr = "!cmake -D EXASIM_NOMPI=OFF -D EXASIM_MPI=ON -D EXASIM_CUDA=OFF ../install";
  end  
end

eval(comstr);  
eval("!cmake --build .");

cd(char(cdir));

