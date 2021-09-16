function cmakecompile(pde, numpde)

if nargin<2
    numpde=1;
end

cdir = pwd(); ii = strfind(cdir, "Exasim");
ExasimPath = cdir(1:(ii+5));
DataPath = ".." + cdir((ii+6):end);
versiondir = ExasimPath + "/" + pde.version;
appdriverdir = versiondir  + "/Kernel/AppDriver/";

copyfile(char(appdriverdir + "opuApp.cpp"), char(cdir + "/app"));
copyfile(char(appdriverdir + "cpuApp.cpp"), char(cdir + "/app"));
copyfile(char(appdriverdir + "gpuApp.cu"), char(cdir + "/app"));
copyfile('app/*.cpp', char(ExasimPath + "/Applications/App"));
copyfile('app/*.cu', char(ExasimPath + "/Applications/App"));

if lower(pde.version) == lower("Version0.1")
    version = " -D EXASIM_VERSION01=ON ";
    pdenum = "";    
elseif lower(pde.version) == lower("Version0.2")
    version = " -D EXASIM_VERSION02=ON ";    
    pdenum = " " + num2str(numpde) + " ";
elseif lower(pde.version) == lower("Version0.3")
    version = " -D EXASIM_VERSION03=ON ";
    pdenum = " " + num2str(numpde) + " ";
end
version = version + "-D EXASIM_APPDIR=" + DataPath + "/app ";

% create exec folder if it does not exist
bindir = "exec";
cd(char(ExasimPath));
if exist(bindir, "file") == 0
    mkdir(char(bindir));    
end
cd(char(ExasimPath + "/" + bindir));

if exist("CMakeCache.txt", "file")
    delete(char("CMakeCache.txt"));
end
if pde.platform == "gpu"
    if pde.buildexec
        eval("!cmake -D EXASIM_APP=ON -D EXASIM_CORES=ON -D EXASIM_EXECUTABLES=ON -D EXASIM_CUDA=ON " + version + "../Installation");
    else
        if exist("libcpuCore.a", "file") && exist("libgpuCore.a", "file") && exist("gpuEXASIM", "file")
            eval("!cmake -D EXASIM_APP=ON -D EXASIM_CUDA=ON " + version + "../Installation");
        else
            eval("!cmake -D EXASIM_APP=ON -D EXASIM_CORES=ON -D EXASIM_EXECUTABLES=ON -D EXASIM_CUDA=ON " + version + "../Installation");
        end
    end
elseif pde.platform == "cpu"
    if pde.buildexec
        eval("!cmake -D EXASIM_APP=ON -D EXASIM_CORES=ON -D EXASIM_EXECUTABLES=ON " + version + "../Installation");
    else
        if exist("libcpuCore.a", "file") && exist("cpuEXASIM", "file")
           eval("!cmake -D EXASIM_APP=ON " + version + "../Installation");
        else
            eval("!cmake -D EXASIM_APP=ON -D EXASIM_CORES=ON -D EXASIM_EXECUTABLES=ON " + version + "../Installation");
        end
    end
end
eval("!cmake --build .");

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

