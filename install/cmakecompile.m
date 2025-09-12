function comstr = cmakecompile(pde,mpiprocs)

if nargin<2
  mpiprocs = pde.mpiprocs;
end

disp("Compile C++ Exasim code using cmake...")

cdir = pwd();
buildpath = char(pde.buildpath);
cd(buildpath);

ii = strfind(pde.buildpath, pde.codename);
n = length(strfind(buildpath(ii:end),'/'));
mystr = "install";
for i = 1:n
  mystr = "../" + mystr;
end

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

if mpiprocs==1 
  if pde.platform == "gpu"
    comstr = "!cmake -D EXASIM_NOMPI=ON -D EXASIM_MPI=OFF -D EXASIM_CUDA=ON " + mystr;    
  elseif pde.platform == "hip"
    comstr = "!cmake -D CMAKE_CXX_COMPILER=hipcc -D EXASIM_NOMPI=ON -D EXASIM_HIP=ON " + mystr;      
  else
    comstr = "!cmake -D EXASIM_NOMPI=ON -D EXASIM_MPI=OFF -D EXASIM_CUDA=OFF " + mystr;
  end
else
  if pde.platform == "gpu"
    comstr = "!cmake -D EXASIM_NOMPI=OFF -D EXASIM_MPI=ON -D EXASIM_CUDA=ON " + mystr;  
  elseif pde.platform == "hip"
    comstr = "!cmake -D CMAKE_CXX_COMPILER=hipcc -D EXASIM_NOMPI=OFF -D EXASIM_MPI=ON -D EXASIM_HIP=ON " + mystr;    
  else
    comstr = "!cmake -D EXASIM_NOMPI=OFF -D EXASIM_MPI=ON -D EXASIM_CUDA=OFF " + mystr;
  end  
end

eval(comstr);  
eval("!cmake --build .");

cd(char(cdir));

