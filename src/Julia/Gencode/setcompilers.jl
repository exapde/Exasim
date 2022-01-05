function setcompilers(app)

cpustatus0 = Sys.which(app.cpucompiler);
cpustatus1 = Sys.which("g++");
cpustatus2 = Sys.which("/usr/bin/g++");
cpustatus3 = Sys.which("/usr/local/bin/g++");
cpustatus4 = Sys.which("/opt/local/bin/g++");

if cpustatus0 != nothing
    print("Using " * app.cpucompiler * " compiler for CPU source code\n");
elseif cpustatus1 != nothing
    app.cpucompiler = "g++";
    print("Using g++ compiler for CPU source code\n");
elseif cpustatus2 != nothing
    app.cpucompiler = "/usr/bin/g++";
    print("Using /usr/bin/g++ compiler for CPU source code\n");
elseif cpustatus3 != nothing
    app.cpucompiler = "/usr/local/bin/g++";
    print("Using /usr/local/bin/g++ compiler for CPU source code\n");
elseif cpustatus4 != nothing
    app.cpucompiler = "/opt/local/bin/g++";
    print("Using /opt/local/bin/g++ compiler for CPU source code\n");
else
    error("Exasim search in /usr/bin, /usr/local/bin, and /opt/local/bin and could not find C++ compiler. Please see the documentation to install it. After installation, please set its path to app.cpucompiler");
end

if app.mpiprocs>1
    mpistatus0 = Sys.which(app.mpicompiler);
    mpistatus1 = Sys.which("mpicxx");
    mpistatus2 = Sys.which("/usr/bin/mpicxx");
    mpistatus3 = Sys.which("/usr/local/bin/mpicxx");
    mpistatus4 = Sys.which("/opt/local/bin/mpicxx");
    mpistatus5 = Sys.which("mpicxx-openmpi-mp");
    mpistatus6 = Sys.which("/usr/bin/mpicxx-openmpi-mp");
    mpistatus7 = Sys.which("/usr/local/bin/mpicxx-openmpi-mp");
    mpistatus8 = Sys.which("/opt/local/bin/mpicxx-openmpi-mp");
    if mpistatus0 != nothing
        print("Using " * app.mpicompiler * " compiler for MPI source code\n");
    elseif mpistatus5 != nothing
        app.mpicompiler = "mpicxx-openmpi-mp";
        app.mpirun = "mpirun-openmpi-mp";
        print("Using mpicxx-openmpi-mp compiler for MPI source code\n");
    elseif mpistatus6 != nothing
        app.mpicompiler = "/usr/bin/mpicxx-openmpi-mp";
        app.mpirun = "/usr/bin/mpirun-openmpi-mp";
        print("Using /usr/bin/mpicxx-openmpi-mp compiler for MPI source code\n");
    elseif mpistatus7 != nothing
        app.mpicompiler = "/usr/local/bin/mpicxx-openmpi-mp";
        app.mpirun = "/usr/local/bin/mpirun-openmpi-mp";
        print("Using /usr/local/bin/mpicxx-openmpi-mp compiler for MPI source code\n");
    elseif mpistatus8 != nothing
        app.mpicompiler = "/opt/local/bin/mpicxx-openmpi-mp";
        app.mpirun = "/opt/local/bin/mpirun-openmpi-mp";
        print("Using opt/local/bin/mpicxx-openmpi-mp compiler for MPI source code\n");
    elseif mpistatus1 != nothing
        app.mpicompiler = "mpicxx";
        app.mpirun = "mpirun";
        print("Using mpicxx compiler for MPI source code\n");
    elseif mpistatus2 != nothing
        app.mpicompiler = "/usr/bin/mpicxx";
        app.mpirun = "/usr/bin/mpirun";
        print("Using /usr/bin/mpicxx compiler for MPI source code\n");
    elseif mpistatus3 != nothing
        app.mpicompiler = "/usr/local/bin/mpicxx";
        app.mpirun = "/usr/local/bin/mpirun";
        print("Using /usr/local/bin/mpicxx compiler for MPI source code\n");
    elseif mpistatus4 != nothing
        app.mpicompiler = "/opt/local/bin/mpicxx";
        app.mpirun = "/opt/local/bin/mpirun";
        print("Using opt/local/bin/mpicxx compiler for MPI source code\n");
    else
        error("Exasim search in /usr/bin, /usr/local/bin, and /opt/local/bin and could not find MPI compiler. Please see the documentation to install it. After installation, please set its path to app.mpicompiler");
    end
    app.mpirun = replace(app.mpicompiler, "mpicxx" => "mpirun");
end

if app.platform == "gpu"
    gpustatus0 = Sys.which(app.gpucompiler);
    gpustatus1 = Sys.which("nvcc");
    gpustatus2 = Sys.which("/usr/bin/nvcc");
    gpustatus3 = Sys.which("/usr/local/bin/nvcc");
    gpustatus4 = Sys.which("/opt/local/bin/nvcc");

    if gpustatus0 != nothing
        print("Using " * app.gpucompiler * " compiler for GPU source code\n");
    elseif gpustatus1 != nothing
        app.gpucompiler = "nvcc";
        print("Using nvcc compiler for GPU source code\n");
    elseif gpustatus2 != nothing
        app.gpucompiler = "/usr/bin/nvcc";
        print("Using /usr/bin/nvcc compiler for GPU source code\n");
    elseif gpustatus3 != nothing
        app.gpucompiler = "/usr/local/bin/nvcc";
        print("Using /usr/local/bin/nvcc compiler for GPU source code\n");
    elseif gpustatus4 != nothing
        app.gpucompiler = "opt/local/bin/nvcc";
        print("Using opt/local/bin/nvcc compiler for GPU source code\n");
    else
        error("Exasim search in /usr/bin, /usr/local/bin, and /opt/local/bin and could not find NVCC compiler. Please see the documentation to install it. After installation, please set its path to app.gpucompiler");
    end
end

mdir = pwd();
ii = findlast(app.codename, mdir);
coredir = mdir[1:ii[end]] * "/lib";

cpulib = 0;
gpulib = 0;
if Sys.isapple()
    if (isfile(coredir * "/Mac/commonCore.a")) && (isfile(coredir * "/Mac/opuCore.a"))
        cpulib = 1;
    end
    if (isfile(coredir * "/Mac/gpuCore.a"))
        gpulib = 1;
    end
elseif Sys.isunix()
    if (isfile(coredir * "/Linux/commonCore.a")) && (isfile(coredir * "/Linux/opuCore.a"))
        cpulib = 1;
    end
    if (isfile(coredir * "/Linux/gpuCore.a"))
        gpulib = 1;
    end
elseif Sys.iswindows()
    if (isfile(coredir * "/Windows/commonCore.a")) && (isfile(coredir * "/Windows/opuCore.a"))
        cpulib = 1;
    end
    if (isfile(coredir * "/Windows/gpuCore.a"))
        gpulib = 1;
    end
end

if cpulib==0
    print("Generating CPU core libraries.");
    genlib(app.cpucompiler, "", coredir);
end

if gpulib==0
    if (length(app.gpucompiler)>0) && (app.platform == "gpu")
        print("Generating GPU core library.");
        genlib("", app.gpucompiler, coredir);
    end
end

return app

end
