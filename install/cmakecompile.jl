function stringcommand(str::String)

ind = findall(" ", str);
n = length(ind);
cmdstr = Array{String,1}(undef,n+1);
cmdstr[1] = str[1:(ind[1][1]-1)];
for i = 2:n
    cmdstr[i] = str[(ind[i-1][1]+1):(ind[i][1]-1)];
end
i = n+1;
cmdstr[i] = str[(ind[i-1][1]+1):end];

return Cmd(cmdstr)

end

function cmakecompile(pde)

display("Compile C++ Exasim code using cmake...")

cdir = pwd(); 
cd(pde.buildpath);

if isfile("cpuEXASIM")
  rm("cpuEXASIM")
end
if isfile("cpumpiEXASIM")
  rm("cpumpiEXASIM")
end
if isfile("gpuEXASIM")
  rm("gpuEXASIM")
end
if isfile("gpumpiEXASIM")
  rm("gpumpiEXASIM")
end

# if isfile("CMakeCache.txt")
#   rm("CMakeCache.txt");
# end

if pde.mpiprocs==1 
  if pde.platform == "gpu"
    comstr = "cmake -D EXASIM_NOMPI=ON -D EXASIM_MPI=OFF -D EXASIM_CUDA=ON ../install";  
  elseif pde.platform == "hip"
    comstr = "cmake -D CMAKE_CXX_COMPILER=hipcc -D EXASIM_NOMPI=ON -D EXASIM_HIP=ON ../install";  
  else
    comstr = "cmake -D EXASIM_NOMPI=ON -D EXASIM_MPI=OFF -D EXASIM_CUDA=OFF ../install";
  end
else
  if pde.platform == "gpu"
    comstr = "cmake -D EXASIM_NOMPI=OFF -D EXASIM_MPI=ON -D EXASIM_CUDA=ON ../install";  
  elseif pde.platform == "hip"
    comstr = "cmake -D CMAKE_CXX_COMPILER=hipcc -D EXASIM_NOMPI=OFF -D EXASIM_MPI=ON -D EXASIM_HIP=ON ../install";  
  else
    comstr = "cmake -D EXASIM_NOMPI=OFF -D EXASIM_MPI=ON -D EXASIM_CUDA=OFF ../install";
  end  
end

run(stringcommand(comstr));
run(stringcommand("cmake --build ."));

cd(cdir);
return comstr;

end

function compilepdemodel(pde)

# Save current working directory
cdir = pwd()

# Construct the path to the build directory
cmp = pde.exasimpath * "/backend/Model/build"

# Change to the build directory
cd(cmp)

# Run cmake and make
run(`cmake ..`)
run(`make`)

# Return to the original working directory
cd(cdir)

return cmp

end


# buildpath = pde.buildpath;
# DataPath = buildpath;
# mpirun = pde.mpirun;
# pdenum = string(numpde) * " ";

# if pde.platform == "cpu"        
#     if pde.mpiprocs==1        
#         str = "./cpuEXASIM " * pdenum * DataPath * "/datain/ " * DataPath * "/dataout/out";        
#     else        
#         str = mpirun * " -np " * string(pde.mpiprocs) * " ./cpumpiEXASIM " * pdenum * DataPath * "/datain/ " * DataPath * "/dataout/out";       
#     end        
# elseif pde.platform == "gpu"
#     if pde.mpiprocs==1        
#         str = "./gpuEXASIM " * pdenum * DataPath * "/datain/ " * DataPath * "/dataout/out";        
#     else        
#         str = mpirun * " -np " * string(pde.mpiprocs) * " ./gpumpiEXASIM " * pdenum * DataPath * "/datain/ " * DataPath * "/dataout/out";       
#     end    
# end

# run(stringcommand(str));

# cd(cdir);

# return str;

