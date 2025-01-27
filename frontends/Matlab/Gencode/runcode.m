function runstr = runcode(pde, numpde, mpiprocs)

if nargin<3
  mpiprocs = pde.mpiprocs;
end

disp("Run C++ Exasim code ...")
cdir = pwd();
cd(char(pde.buildpath));

DataPath = pde.buildpath;

mystr = " " + num2str(numpde) + " ";
if numpde>100 % two-domain problems
  numpde = 2;
end
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
    if mpiprocs==1        
        exec = "./cpuEXASIM ";        
        runstr = "!" + exec + mystr;
    else        
        exec = " ./cpumpiEXASIM ";
        runstr = "!" + mpirun + " -np " + string(mpiprocs) + exec + mystr;
    end        
elseif pde.platform == "gpu"
    if mpiprocs==1        
        exec = "./gpuEXASIM ";     
        runstr = "!" + exec + mystr;
    else        
        exec = " ./gpumpiEXASIM ";
        runstr = "!" + mpirun + " -np " + string(mpiprocs) + exec + mystr;
    end
end

tic
eval(char(runstr));
toc

cd(char(cdir));

end

