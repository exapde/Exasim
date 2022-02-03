function compilecode(app)

print("compile code...\n");

cd("app");

codename = app.codename;
version = app.version;
appname = app.appname;
cpucompiler = app.cpucompiler;
mpicompiler = app.mpicompiler;
gpucompiler = app.gpucompiler;
enzyme = app.enzyme;
cpuflags = app.cpuflags;
gpuflags = app.gpuflags;
cpuappflags = app.cpuappflags
gpuappflags = app.gpuappflags

# current directory
cdir = pwd();
ii = findlast(codename, cdir);
#codedir = cdir[1:ii[end]];

tmp = findall("/", cdir[(ii[end]+1):end]);
up = length(tmp);

codedir = "";
for i = 1:up
    codedir = codedir * "../";
end

if Sys.isapple()
    coredir = codedir * "lib/Mac/";
elseif Sys.isunix()
    coredir = codedir * "lib/Linux/";
elseif Sys.windows()
    coredir = codedir * "lib/Windows/";
end

versiondir = codedir  * version;
appdriverdir = versiondir * "/Kernel/AppDriver/";
maindir = versiondir * "/Kernel/Main/";

cp(appdriverdir * "opuApp.cpp", cdir * "/opuApp.cpp", force=true);
cp(appdriverdir * "cpuApp.cpp", cdir * "/cpuApp.cpp", force=true);
cp(appdriverdir * "gpuApp.cu", cdir * "/gpuApp.cu", force=true);

compilerstr = Array{String, 1}(undef, 12);

if length(cpucompiler)>0
    if (length(enzyme)>0)   
        compilerstr[1] = cpucompiler * " -D _ENZYME -fPIC -O3 -c opuApp.cpp" * " -stdlib=libc++ -Xclang -load -Xclang " * coredir * enzyme;
    else
        compilerstr[1] = cpucompiler * " -fPIC -O3 -c opuApp.cpp";
    end    
    if length(cpuappflags) > 0
        compilerstr[1] = compilerstr[1] * " " * cpuappflags
    end
    compilerstr[2] = "ar -rvs opuApp.a opuApp.o";
else
    compilerstr[1] = "";
    compilerstr[2] = "";
end

if length(gpucompiler)>0
    if (length(enzyme))>0
        compilerstr[3] = gpucompiler * " -D _FORCE_INLINES -O3 -c -fPIC gpuApp.cu";
        compilerstr[3] = compilerstr[3] * " -D _ENZYME -std=c++11 -stdlib=libc++ -Xclang -load -Xclang " * coredir * enzyme;
    else
        compilerstr[3] = gpucompiler * " -D_FORCE_INLINES -O3 -c --compiler-options '-fPIC' gpuApp.cu";
    end
    # compilerstr[3] = compilerstr[3] * " --cuda-gpu-arch=sm_60 -lcudart -L/usr/local/cuda-9.0/lib64 -stdlib=libc++ -std=c++11 --cuda-path=/usr/local/cuda-9.0/"
    if length(gpuappflags) > 0
        compilerstr[3] = compilerstr[3] * " " * gpuappflags
    end
    compilerstr[4] = "ar -rvs gpuApp.a gpuApp.o";
else
    compilerstr[3] = "";
    compilerstr[4] = "";
end

if (length(cpuflags)>0) && (length(cpucompiler)>0)
    #str1 = cpucompiler * " -std=c++11 " * maindir * "main.cpp " * "-o serial" * appname * " ";
    str2 = coredir * "commonCore.a " * coredir * "opuCore.a " * "opuApp.a ";
    str3 = cpuflags;
    #compilerstr[5] = str1 * str2 * str3;
    if (length(enzyme)>0)    
        str1 = cpucompiler  * " -std=c++11 -stdlib=libc++ -D _ENZYME " * maindir * "main.cpp " * "-o serial" * appname * " ";
        compilerstr[5] = str1 * str2 * str3 * " -Xclang -load -Xclang " * coredir * enzyme;
    else
        str1 = cpucompiler * " -std=c++11 " * maindir * "main.cpp " * "-o serial" * appname * " ";
        compilerstr[5] = str1 * str2 * str3;
    end
end

if (length(cpuflags)>0) && (length(mpicompiler)>0)
    if (length(enzyme)>0)
        str1 = mpicompiler * " -std=c++11 -D _ENZYME -D _MPI " * maindir * "main.cpp " * "-o mpi" * appname * " ";
    else
        str1 = mpicompiler * " -std=c++11 -D _MPI " * maindir * "main.cpp " * "-o mpi" * appname * " ";
    end
    str2 = coredir * "commonCore.a " * coredir * "opuCore.a " * "opuApp.a ";
    str3 = cpuflags;
    compilerstr[6] = str1 * str2 * str3;
end

if (length(cpuflags)>0) && (length(cpucompiler)>0) && (length(gpucompiler)>0) && (length(gpuflags)>0)
    if (length(enzyme)>0)
        str1 = cpucompiler * " -std=c++11 -D _ENZYME -D _CUDA " * maindir * "main.cpp " * "-o gpu" * appname * " ";
        str1 = str1 * " -Xclang -load -Xclang " * coredir * enzyme * " ";
    else
        str1 = cpucompiler * " -std=c++11 -D _CUDA " * maindir * "main.cpp " * "-o gpu" * appname * " ";
    end
    str2 = coredir * "commonCore.a " * coredir * "gpuCore.a " * coredir * "opuCore.a opuApp.a gpuApp.a ";
    str3 = cpuflags * " " * gpuflags;
    compilerstr[7] = str1 * str2 * str3;
end

if (length(cpuflags)>0) && (length(mpicompiler)>0) && (length(gpucompiler)>0) && (length(gpuflags)>0)
    str1 = mpicompiler * " -std=c++11  -D _MPI -D _CUDA " * maindir * "main.cpp " * "-o gpumpi" * appname * " ";
    str2 = coredir * "commonCore.a " * coredir * "gpuCore.a " * coredir * "opuCore.a opuApp.a gpuApp.a ";
    str3 = cpuflags * " " * gpuflags;
    compilerstr[8] = str1 * str2 * str3;
end

if length(cpucompiler)>0
    compilerstr[9] = cpucompiler * " -fPIC -O3 -c cpuApp.cpp -fopenmp";
    compilerstr[10] = "ar -rvs cpuApp.a cpuApp.o";
else
    compilerstr[9] = "";
    compilerstr[10] = "";
end

if (length(cpuflags)>0) && (length(cpucompiler)>0)
    str1 = cpucompiler * " -std=c++11 " * maindir * "main.cpp" * "-o openmp" * appname * " ";
    str2 = coredir * "commonCore.a " * coredir * "cpuCore.a cpuApp.a "
    str3 = "-fopenmp " * cpuflags;
    compilerstr[11] = str1 * str2 * str3;
end

if (length(cpuflags)>0) && (length(mpicompiler)>0)
    str1 = mpicompiler * " -std=c++11 -D _MPI " * maindir * "main.cpp" * "-o openmpmpi" * appname * " ";
    str2 = coredir * "commonCore.a " * coredir * "cpuCore.a cpuApp.a ";
    str3 = "-fopenmp " * cpuflags;
    compilerstr[12] = str1 * str2 * str3;
end

if app.platform == "cpu"
    run(string2cmd(compilerstr[1]));
    run(string2cmd(compilerstr[2]));
    if app.mpiprocs==1
        run(string2cmd(compilerstr[5]));
    else
        run(string2cmd(compilerstr[6]));
    end
elseif app.platform == "gpu"
    println(compilerstr[1])
    println(compilerstr[2])
    println(compilerstr[3])
    println(compilerstr[4])

    run(string2cmd(compilerstr[1]));
    run(string2cmd(compilerstr[2]));
    run(string2cmd(compilerstr[3]));
    run(string2cmd(compilerstr[4]));
    if app.mpiprocs==1
        println(compilerstr[7])
        run(string2cmd(compilerstr[7]));
    else
        println(compilerstr[8])
        run(string2cmd(compilerstr[8]));
    end
end

cd("..");

return compilerstr;

end
