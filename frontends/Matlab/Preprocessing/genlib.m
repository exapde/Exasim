function genlib(cpucompiler, gpucompiler, coredir, cpulibflags, gpulibflags)
% Examples: genlib("g++"); 
%           genlib("g++","nvcc"); 
if nargin < 4
    cpulibflags = "";
end
if nargin < 5 
    gpulibflags = "";
end

mydir = pwd;
cd(coredir);

if ~isempty(cpucompiler)
    str = "!" + cpucompiler + " -fPIC -O3 -c commonCore.cpp";
    str = str + " " + cpulibflags;
    eval(char(str));        
    !ar rvs commonCore.a commonCore.o

    str = "!" + cpucompiler + " -fPIC -O3 -c opuCore.cpp";
    str = str + " " + cpulibflags;
    eval(char(str));    
    !ar rvs opuCore.a opuCore.o

    str = "!" + cpucompiler + " -fPIC -O3 -c cpuCore.cpp -fopenmp";
    str = str + " " + cpulibflags;
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
    str = str + " " + gpulibflags;
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







