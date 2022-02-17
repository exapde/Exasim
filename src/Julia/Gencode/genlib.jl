function genlib(cpucompiler::String, gpucompiler::String, coredir::String, cpulibflags::String="", gpulibflags::String="")
# Examples: genlib("g++","");
#           genlib("g++","nvcc");

mydir = pwd();
cd(coredir);

if length(cpucompiler)>0
    str = cpucompiler * " -fPIC -O3 -c commonCore.cpp"
    if length(cpulibflags) > 0
        str = str * " " * cpulibflags
    end
    run(string2cmd(str));

    str = `ar rvs commonCore.a commonCore.o`;
    run(str);

    str = cpucompiler * " -fPIC -O3 -c opuCore.cpp"
    if length(cpulibflags)
        str = str * " " * cpulibflags
    end
    run(string2cmd(str));

    str = `ar rvs opuCore.a opuCore.o`;
    run(str);

    # str = `$cpucompiler -fPIC -O3 -c cpuCore.cpp -fopenmp`;
    # run(str);
    # str = `ar rvs cpuCore.a cpuCore.o`;
    # run(str);
end

if length(gpucompiler)>0
    str = gpucompiler * "-D_FORCE_INLINES -O3 -c --compiler-options '-fPIC' gpuCore.cu";
    if length(gpulibflags) > 0
        str = str * " " * gpulibflags
    end
    run(string2cmd(str));
    str = `ar rvs gpuCore.a gpuCore.o `;
    run(str);
end

if Sys.isapple()
#    mv("*.a",pwd() * "/Mac");
    if !isdir("Mac")
        mkdir("Mac")
    end
    if isfile("commonCore.a")
        run(`mv commonCore.a Mac`);
    end
    if isfile("opuCore.a")
        run(`mv opuCore.a Mac`);
    end
    if isfile("cpuCore.a")
        run(`mv cpuCore.a Mac`);
    end
    if isfile("gpuCore.a")
        run(`mv gpuCore.a Mac`);
    end
elseif Sys.isunix()
    if !isdir("Linux")
        mkdir("Linux")
    end
    if isfile("commonCore.a")
        run(`mv commonCore.a Linux`);
    end
    if isfile("opuCore.a")
        run(`mv opuCore.a Linux`);
    end
    if isfile("cpuCore.a")
        run(`mv cpuCore.a Linux`);
    end
    if isfile("gpuCore.a")
        run(`mv gpuCore.a Linux`);
    end
elseif Sys.iswindows()
    if !isdir("Windows")
        mkdir("Windows")
    end
    if isfile("commonCore.a")
        run(`mv commonCore.a Windows`);
    end
    if isfile("opuCore.a")
        run(`mv opuCore.a Windows`);
    end
    if isfile("cpuCore.a")
        run(`mv cpuCore.a Windows`);
    end
    if isfile("gpuCore.a")
        run(`mv gpuCore.a Windows`);
    end
end

rm("commonCore.a", force=true);
rm("opuCore.a", force=true);
rm("cpuCore.a", force=true);
rm("gpuCore.a", force=true);
rm("commonCore.o", force=true);
rm("opuCore.o", force=true);
rm("cpuCore.o", force=true);
rm("gpuCore.o", force=true);

cd(mydir);

end
