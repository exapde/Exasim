function compilerstr = compilecode(app)

disp("compile code...");

cd(char("app"));
codename = app.codename;
version = app.version;
appname = app.appname;
cpucompiler = app.cpucompiler;
mpicompiler = app.mpicompiler;
gpucompiler = app.gpucompiler;
cpuflags = app.cpuflags;
gpuflags = app.gpuflags;

% current directory
cdir = pwd();
ii = strfind(cdir, codename);

tmp = strfind(cdir(ii(end):end), "/");
up = length(tmp);
codedir = "";
for i = 1:up
    codedir = codedir + "../";
end

if ismac
    coredir = codedir + "Core/Mac/";
elseif isunix
    coredir = codedir + "Core/Linux/";
elseif ispc
    coredir = codedir + "Core/Windows/";
end

versiondir = codedir  + version;
appdriverdir = versiondir + "/Kernel/AppDriver/";
maindir = versiondir + "/Kernel/Main/";

copyfile(char(appdriverdir + "opuApp.cpp"), char("opuApp.cpp"));
copyfile(char(appdriverdir + "cpuApp.cpp"), char("cpuApp.cpp"));
copyfile(char(appdriverdir + "gpuApp.cu"), char("gpuApp.cu"));

compilerstr = cell(12,1);
for i = 1:12
    compilerstr{i} = "";
end

if ~isempty(cpucompiler)
    compilerstr{1} = cpucompiler + " -fPIC -O3 -c opuApp.cpp";
    compilerstr{2} = "ar -rvs opuApp.a opuApp.o";
end

if ~isempty(gpucompiler)
    compilerstr{3} = gpucompiler + " -D_FORCE_INLINES -O3 -c --compiler-options '-fPIC' -w gpuApp.cu";
    compilerstr{4} = "ar -rvs gpuApp.a gpuApp.o";
end

if (~isempty(cpuflags)>0) && (~isempty(cpucompiler)>0)
    str1 = cpucompiler + " -std=c++11 " + maindir + "main.cpp " + "-o serial" + appname + " ";
    str2 = coredir + "commonCore.a " + coredir + "opuCore.a " + "opuApp.a ";
    str3 = cpuflags;
    compilerstr{5} = str1 + str2 + str3;
end

if (~isempty(cpuflags)) && (~isempty(mpicompiler)>0)
    str1 = mpicompiler + " -std=c++11 -D _MPI " + maindir + "main.cpp " + "-o mpi" + appname + " ";
    str2 = coredir + "commonCore.a " + coredir + "opuCore.a " + "opuApp.a ";
    str3 = cpuflags;
    compilerstr{6} = str1 + str2 + str3;
end

if (~isempty(cpuflags)>0) && (~isempty(cpucompiler)>0) && (~isempty(gpucompiler)>0) && (~isempty(gpuflags)>0)
    str1 = cpucompiler + " -std=c++11 -D _CUDA " + maindir + "main.cpp " + "-o gpu" + appname + " ";
    str2 = coredir + "commonCore.a " + coredir + "gpuCore.a " + coredir + "opuCore.a opuApp.a gpuApp.a ";
    str3 = cpuflags + " " + gpuflags;
    compilerstr{7} = str1 + str2 + str3;
end

if (~isempty(cpuflags)) && (~isempty(mpicompiler)>0) && (~isempty(gpucompiler)>0) && (~isempty(gpuflags)>0)
    str1 = mpicompiler + " -std=c++11  -D _MPI -D _CUDA " + maindir + "main.cpp " + "-o gpumpi" + appname + " ";
    str2 = coredir + "commonCore.a " + coredir + "gpuCore.a " + coredir + "opuCore.a opuApp.a gpuApp.a ";
    str3 = cpuflags + " " + gpuflags;
    compilerstr{8} = str1 + str2 + str3;
end

if ~isempty(cpucompiler)
    compilerstr{9} = cpucompiler + " -fPIC -O3 -c cpuApp.cpp -fopenmp";
    compilerstr{10} = "ar -rvs cpuApp.a cpuApp.o";
end

if (~isempty(cpuflags)>0) && (~isempty(cpucompiler))
    str1 = cpucompiler + " -std=c++11 " + maindir + "main.cpp" + "-o openmp" + appname + " ";
    str2 = coredir + "commonCore.a " + coredir + "cpuCore.a cpuApp.a ";
    str3 = "-fopenmp " + cpuflags;
    compilerstr{11} = str1 + str2 + str3;
end

if (~isempty(cpuflags)) && (~isempty(mpicompiler))
    str1 = mpicompiler + " -std=c++11 -D _MPI " + maindir + "main.cpp" + "-o openmpmpi" + appname + " ";
    str2 = coredir + "commonCore.a " + coredir + "cpuCore.a cpuApp.a ";
    str3 = "-fopenmp " + cpuflags;
    compilerstr{12} = str1 + str2 + str3;
end

delete *.a;

if app.platform == "cpu"
    eval(char("!" + compilerstr{1}));
    eval(char("!" + compilerstr{2}));    
    if app.mpiprocs==1
        eval(char("!" + compilerstr{5}));
    else
        eval(char("!" + compilerstr{6}));
    end
elseif app.platform == "gpu"
    eval(char("!" + compilerstr{1}));
    eval(char("!" + compilerstr{2}));        
    eval(char("!" + compilerstr{3}));
    eval(char("!" + compilerstr{4}));
    if app.mpiprocs==1
        eval(char("!" + compilerstr{7}));
    else
        eval(char("!" + compilerstr{8}));
    end
end

delete *.o;

cd(char(".."));

end

