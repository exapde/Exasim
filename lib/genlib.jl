function genlib(cpucompiler::String, gpucompiler::String)
# Examples: genlib("g++","");
#           genlib("g++","nvcc");

if length(cpucompiler)==0
    error("cpucompiler is empty");
end

str = `$cpucompiler -fPIC -O3 -c commonCore.cpp`;
run(str);
str = `ar rvs commonCore.a commonCore.o`;
run(str);

str = `$cpucompiler -fPIC -O3 -c opuCore.cpp`;
run(str);
str = `ar rvs opuCore.a opuCore.o`;
run(str);

# str = `$cpucompiler -fPIC -O3 -c cpuCore.cpp -fopenmp`;
# run(str);
# str = `ar rvs cpuCore.a cpuCore.o`;
# run(str);

if length(gpucompiler)>0
    str = `$gpucompiler -D_FORCE_INLINES -O3 -c --compiler-options '-fPIC' gpuCore.cu`;
    run(str);
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

end

# delete *.o;
# delete *.a;
