function genlib(cpucompiler, gpucompiler, coredir)
% Examples: genlib("g++"); 
%           genlib("g++","nvcc"); 

mydir = pwd;
cd(coredir);

if ~isempty(cpucompiler)
    str = "!" + cpucompiler + " -fPIC -O3 -c commonCore.cpp";
    eval(char(str));        
    !ar rvs commonCore.a commonCore.o

    str = "!" + cpucompiler + " -fPIC -O3 -c opuCore.cpp";
    eval(char(str));    
    !ar rvs opuCore.a opuCore.o

    str = "!" + cpucompiler + " -fPIC -O3 -c cpuCore.cpp -fopenmp";
    eval(char(str));    
    !ar rvs cpuCore.a cpuCore.o   
    
    if ismac    
        mkdir("Mac");
        copyfile('*.a','Mac');    
    elseif isunix
        mkdir("Linux");
        copyfile('*.a','Linux');
    elseif ispc
        mkdir("Windows");
        copyfile('*.a','Windows');
    end
   
    delete *.o;
    delete *.a;    
end

if ~isempty(gpucompiler)
    str = "!" + gpucompiler + " -D_FORCE_INLINES -O3 -c --compiler-options '-fPIC' gpuCore.cu";
    eval(char(str));    
    !ar -rvs gpuCore.a gpuCore.o  
    
    if ismac    
        copyfile('*.a','Mac');    
    elseif isunix
        copyfile('*.a','Linux');
    elseif ispc
        copyfile('*.a','Windows');
    end
   
    delete *.o;
    delete *.a;        
end

cd(mydir);







