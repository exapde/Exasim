function checkcompilers(app)

cpustatus = Sys.which(app.cpucompiler);
if cpustatus != nothing
    display("Using " * app.cpucompiler * " compiler for CPU source code");
else
    display(app.cpucompiler * " compiler is not available on your system.");
    cpustatus = Sys.which("g++");
    if cpustatus != nothing
        display("However, g++ compiler is available on your system. It will be used for compiling CPU source code.");
        app.cpucompiler = "g++";
    else
        error("C++ compiler is not available on your system.");
    end
end

if app.mpiprocs>1
    mpistatus = Sys.which(app.mpicompiler);
    if mpistatus != nothing
        display("Using " * app.mpicompiler * " compiler for MPI source code");
    else
        display(app.mpicompiler * " compiler is not available on your system.");
        mpistatus = Sys.which("mpicxx");
        if mpistatus != nothing
            display("However, mpicxx compiler is available on your system. It will be used for compiling MPI source code.");
            app.mpicompiler = "mpicxx";
            app.mpirun = "mpirun";
        else
            error("MPI compiler is not available on your system.");
        end
    end
end

if app.platform == "gpu"
    gpustatus = Sys.which(app.gpucompiler);
    if gpustatus != nothing
        display("Using " * app.gpucompiler * " compiler for GPU source code");
    else
        display(app.gpucompiler * " compiler is not available on your system.");
        gpustatus = Sys.which("nvcc");
        if gpustatus != nothing
            display("However, nvcc compiler is available on your system. It will be used for compiling GPU source code.");
            app.gpucompiler = "nvcc";
        else
            error("GPU compiler is not available on your system.");
        end
    end
end

mdir = pwd();
ii = findlast(app.codename, mdir);
coredir = mdir[1:ii[end]] * "/Core";

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
    display("Generating CPU core libraries.");
    genlib(app.cpucompiler, "", coredir);
end

if gpulib==0
    if (length(app.gpucompiler)>0) && (app.platform == "gpu")
        display("Generating GPU core library.");
        genlib("", app.gpucompiler, coredir);
    end
end

return cpulib, gpulib

end
