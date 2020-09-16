function [app,cpulib,gpulib] = setcompilers(app)

% app.gmsh = findexec(app.gmsh, app.version);
% app.paraview = findexec(app.paraview, app.version);
% app.metis = findexec(app.metis, app.version);

[cpustatus0,~] = system("which " + app.cpucompiler);
[cpustatus1,~] = system("which g++");
[cpustatus2,~] = system("which /usr/bin/g++");
[cpustatus3,~] = system("which /usr/local/bin/g++");
[cpustatus4,~] = system("which /opt/local/bin/g++");        

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
    error("Exasim search in /usr/bin, /usr/local/bin, and /opt/local/bin and could not find C++ compiler. Please see the documentation to install it. After installation, please set its path to app.cpucompiler"); 
end

if app.mpiprocs>1
    [mpistatus0,~] = system("which "+ app.mpicompiler);
    [mpistatus1,~] = system("which mpicxx");
    [mpistatus2,~] = system("which /usr/bin/mpicxx");
    [mpistatus3,~] = system("which /usr/local/bin/mpicxx");
    [mpistatus4,~] = system("which /opt/local/bin/mpicxx");        
    [mpistatus5,~] = system("which mpicxx-openmpi-mp");
    [mpistatus6,~] = system("which /usr/bin/mpicxx-openmpi-mp");
    [mpistatus7,~] = system("which /usr/local/bin/mpicxx-openmpi-mp");
    [mpistatus8,~] = system("which /opt/local/bin/mpicxx-openmpi-mp");        
    
    if mpistatus0==0
        disp("Using " + app.mpicompiler + " compiler for MPI source code");
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
    else            
        error("Exasim search in /usr/bin, /usr/local/bin, and /opt/local/bin and could not find MPI compiler. Please see the documentation to install it. After installation, please set its path to app.mpicompiler");         
    end        
    app.mpirun = strrep(app.mpicompiler, "mpicxx", "mpirun");
end

if app.platform == "gpu"
    [gpustatus0,~] = system("which " + app.gpucompiler);
    [gpustatus1,~] = system("which nvcc");
    [gpustatus2,~] = system("which /usr/bin/nvcc");
    [gpustatus3,~] = system("which /usr/local/bin/nvcc");
    [gpustatus4,~] = system("which /opt/local/bin/nvcc");        

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
        error("Exasim search in /usr/bin, /usr/local/bin, and /opt/local/bin and could not find NVCC compiler. Please see the documentation to install it. After installation, please set its path to app.gpucompiler");         
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

