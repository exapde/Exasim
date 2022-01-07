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

function cmakecompile(pde, numpde=1)

cdir = pwd(); ii = findlast("Exasim", cdir);
ExasimPath = cdir[1:ii[end]];
DataPath = ".." * cdir[(ii[end]+1):end];
versiondir = ExasimPath * "/" * pde.version;
appdriverdir = versiondir  * "/Kernel/AppDriver/";

cp(appdriverdir * "opuApp.cpp", cdir * "/app/opuApp.cpp", force=true);
cp(appdriverdir * "cpuApp.cpp", cdir * "/app/cpuApp.cpp", force=true);
cp(appdriverdir * "gpuApp.cu",  cdir * "/app/gpuApp.cu", force=true);
#cp("app/*.cpp", ExasimPath * "/Applications/App", force=true);
#cp("app/*.cu",  ExasimPath * "/Applications/App", force=true);
#run(stringcommand("cp -f " * "app/*.cpp " * ExasimPath * "/Applications/App"));
#run(stringcommand("cp -f " * "app/*.cu " * ExasimPath * "/Applications/App"));

if lowercase(pde.version) == lowercase("Version0.1")
    version = " -D EXASIM_VERSION01=ON ";
    pdenum = "";    
elseif lowercase(pde.version) == lowercase("Version0.2")
    version = " -D EXASIM_VERSION02=ON ";    
    pdenum = string(numpde) * " ";
elseif lowercase(pde.version) == lowercase("Version0.3")
    version = " -D EXASIM_VERSION03=ON ";
    pdenum = string(numpde) * " ";
end
version = version * "-D EXASIM_APPDIR=" * DataPath * "/app ";

# create exec folder if it does not exist
bindir = "exec";
cd(ExasimPath);
if !isdir(bindir)    
    mkdir(bindir);    
end
cd(ExasimPath * "/" * bindir);

if isfile("CMakeCache.txt")
    rm("CMakeCache.txt");
end
if pde.platform == "gpu"
    if pde.buildexec==1
        run(stringcommand("cmake -D EXASIM_APP=ON -D EXASIM_CORES=ON -D EXASIM_EXECUTABLES=ON -D EXASIM_CUDA=ON " * version * "../Installation"));
    else
        if exist("libcpuCore.a", "file") && exist("libgpuCore.a", "file") && exist("gpuEXASIM", "file")
            run(stringcommand("cmake -D EXASIM_APP=ON -D EXASIM_CUDA=ON " * version * "../Installation"));
        else
            run(stringcommand("cmake -D EXASIM_APP=ON -D EXASIM_CORES=ON -D EXASIM_EXECUTABLES=ON -D EXASIM_CUDA=ON " * version * "../Installation"));
        end
    end
elseif pde.platform == "cpu"
    if pde.buildexec==1
        run(stringcommand("cmake -D EXASIM_APP=ON -D EXASIM_CORES=ON -D EXASIM_EXECUTABLES=ON " * version * "../Installation"));
    else
        if exist("libcpuCore.a", "file") && exist("cpuEXASIM", "file")
           run(stringcommand("cmake -D EXASIM_APP=ON " * version * "../Installation"));
        else
           run(stringcommand("cmake -D EXASIM_APP=ON -D EXASIM_CORES=ON -D EXASIM_EXECUTABLES=ON " * version * "../Installation"))s;
        end
    end
end
run(stringcommand("cmake --build ."));

mpirun = pde.mpirun;
if pde.platform == "cpu"    
    if pde.mpiprocs==1
        str = "./cpuEXASIM " * pdenum * DataPath * "/datain/ " * DataPath * "/dataout/out";        
    else
        str = mpirun * " -np " * string(pde.mpiprocs) * " ./cpumpiEXASIM " * pdenum * DataPath * "/datain/ " * DataPath * "/dataout/out";       
    end    
    run(stringcommand(str));
elseif pde.platform == "gpu"
    if pde.mpiprocs==1
        str = "./gpuEXASIM " * pdenum * DataPath * "/datain/ " * DataPath * "/dataout/out";        
    else
        str = mpirun * " -np " * string(pde.mpiprocs) * " ./app/gpumpiEXASIM " * pdenum * DataPath * "/datain/ " * DataPath * "/dataout/out";       
    end    
    run(stringcommand(str));
end

cd(cdir);

return str;

end

