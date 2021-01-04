function [app,cpulib,gpulib] = checkcompilers(app)

[cpustatus0,~] = system(app.cpucompiler + " -v");
[cpustatus1,~] = system("g++ -v");
[cpustatus2,~] = system("/usr/bin/g++ -v");
[cpustatus3,~] = system("/usr/local/bin/g++ -v");
[cpustatus4,~] = system("/opt/local/bin/g++ -v");        

if cpustatus0==0
    disp("Using " + app.cpucompiler + " compiler for CPU source code");
elseif cpustatus1==0        
    app.cpucompiler = "g++";
    disp("Using g++ compiler for CPU source code");
elseif cpustatus2==0        
    app.cpucompiler = "/usr/bin/g++";    
    disp("Using /usr/bin/g++ compiler for CPU source code");
elseif cpustatus3==0        
    app.cpucompiler = "/usr/local/bin/g++";
    disp("Using /usr/local/bin/g++ compiler for CPU source code");
elseif cpustatus4==0        
    app.cpucompiler = "/opt/local/bin/g++";    
    disp("Using /opt/local/bin/g++ compiler for CPU source code");
else            
    error("C++ compiler is not found on your system. Please install it. If a C++ compiler is available on your system, please set its path to app.cpucompiler"); 
end

if app.mpiprocs>1
    [mpistatus0,~] = system(app.mpicompiler + " -v");
    [mpistatus1,~] = system("mpicxx -v");
    [mpistatus2,~] = system("/usr/bin/mpicxx -v");
    [mpistatus3,~] = system("/usr/local/bin/mpicxx -v");
    [mpistatus4,~] = system("/opt/local/bin/mpicxx -v");        
    [mpistatus5,~] = system("mpicxx-openmpi-mp -v");
    [mpistatus6,~] = system("/usr/bin/mpicxx-openmpi-mp -v");
    [mpistatus7,~] = system("/usr/local/bin/mpicxx-openmpi-mp -v");
    [mpistatus8,~] = system("/opt/local/bin/mpicxx-openmpi-mp -v");        
    
    if mpistatus0==0
        disp("Using " + app.mpicompiler + " compiler for MPI source code");
    elseif mpistatus5==0        
        app.mpicompiler = "mpicxx-openmpi-mp";
        app.mpirun = "mpirun-openmpi-mp";
        disp("Using mpicxx-openmpi-mp compiler for MPI source code");
    elseif mpistatus6==0        
        app.mpicompiler = "/usr/bin/mpicxx-openmpi-mp";    
        app.mpirun = "/usr/bin/mpirun-openmpi-mp";
        disp("Using /usr/bin/mpicxx-openmpi-mp compiler for MPI source code");
    elseif mpistatus7==0        
        app.mpicompiler = "/usr/local/bin/mpicxx-openmpi-mp";    
        app.mpirun = "/usr/local/bin/mpirun-openmpi-mp";
        disp("Using /usr/local/bin/mpicxx-openmpi-mp compiler for MPI source code");
    elseif mpistatus8==0        
        app.mpicompiler = "/opt/local/bin/mpicxx-openmpi-mp";    
        app.mpirun = "/opt/local/bin/mpirun-openmpi-mp";                
        disp("Using opt/local/bin/mpicxx-openmpi-mp compiler for MPI source code");                
    elseif mpistatus1==0        
        app.mpicompiler = "mpicxx";
        app.mpirun = "mpirun";
        disp("Using mpicxx compiler for MPI source code");
    elseif mpistatus2==0        
        app.mpicompiler = "/usr/bin/mpicxx";    
        app.mpirun = "/usr/bin/mpirun";
        disp("Using /usr/bin/mpicxx compiler for MPI source code");
    elseif mpistatus3==0        
        app.mpicompiler = "/usr/local/bin/mpicxx";    
        app.mpirun = "/usr/local/bin/mpirun";
        disp("Using /usr/local/bin/mpicxx compiler for MPI source code");
    elseif mpistatus4==0        
        app.mpicompiler = "/opt/local/bin/mpicxx";    
        app.mpirun = "/opt/local/bin/mpirun";                
        disp("Using opt/local/bin/mpicxx compiler for MPI source code");
    else            
        error("MPI compiler is not found on your system. Please install it. If a MPI compiler is available on your system, please set its path to app.mpicompiler"); 
    end        
end

if app.platform == "gpu"
    [gpustatus0,~] = system(app.gpucompiler + " -v");
    [gpustatus1,~] = system("nvcc -v");
    [gpustatus2,~] = system("/usr/bin/nvcc -v");
    [gpustatus3,~] = system("/usr/local/bin/nvcc -v");
    [gpustatus4,~] = system("/opt/local/bin/nvcc -v");        

    if gpustatus0==0
        disp("Using " + app.gpucompiler + " compiler for GPU source code");
    elseif gpustatus1==0        
        app.gpucompiler = "nvcc";
        disp("Using nvcc compiler for GPU source code");
    elseif gpustatus2==0        
        app.gpucompiler = "/usr/bin/nvcc";    
        disp("Using /usr/bin/nvcc compiler for GPU source code");
    elseif gpustatus3==0        
        app.gpucompiler = "/usr/local/bin/nvcc";
        disp("Using /usr/local/bin/nvcc compiler for GPU source code");
    elseif gpustatus4==0        
        app.gpucompiler = "opt/local/bin/nvcc";    
        disp("Using opt/local/bin/nvcc compiler for GPU source code");
    else            
        error("GPU compiler is not found on your system. Please install it. If a GPU compiler is available on your system, please set its path to app.gpucompiler"); 
    end    
end

mdir = pwd(); ii = strfind(mdir, "Exasim");
coredir = mdir(1:(ii+5)) + "/Core";

cpulib = 0;
gpulib = 0;
if ismac    
    if (exist(coredir + "/Mac/commonCore.a", 'file') == 2) && (exist(coredir + "/Mac/opuCore.a", 'file') == 2)
        cpulib = 1;
    end
    if (exist(coredir + "/Mac/gpuCore.a", 'file') == 2)
        gpulib = 1;
    end
elseif isunix
    if (exist(coredir + "/Linux/commonCore.a", 'file') == 2) && (exist(coredir + "/Linux/opuCore.a", 'file') == 2)
        cpulib = 1;
    end
    if (exist(coredir + "/Linux/gpuCore.a", 'file') == 2)
        gpulib = 1;
    end
elseif ispc
    if (exist(coredir + "/Windows/commonCore.a", 'file') == 2) && (exist(coredir + "/Windows/opuCore.a", 'file') == 2)
        cpulib = 1;
    end
    if (exist(coredir + "/Windows/gpuCore.a", 'file') == 2)
        gpulib = 1;
    end
end

if cpulib==0    
    disp("Generating CPU core libraries.");
    genlib(app.cpucompiler, [], coredir);
end

if gpulib==0        
    if (~isempty(app.gpucompiler)) && (app.platform == "gpu")
        disp("Generating GPU core library.");
        genlib([], app.gpucompiler, coredir);
    end
end

end

