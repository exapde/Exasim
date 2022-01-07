function compilerstr = compilecppcode(version,appname,cpucompiler,mpicompiler,gpucompiler,cpuflags,gpuflags,blaslapack,ompblaslapack,deletesourcefiles)

if isempty(blaslapack)
    % using default blas and lapack library
    blaslapack = '-lblas -llapack';
end

if isempty(cpuflags)
    % using default blas and lapack library
    cpuflags = '';
end

if ispc
    sslash = '\';
elseif isunix
    sslash = '/';
end
ps=strcat(pwd,sslash);
is=find(ps==sslash);
up=0;
while ~(strcmp(ps(is(end-up-1)+1:is(end-up)-1),'exasim') || strcmp(ps(is(end-up-1)+1:is(end-up)-1),'EXASIM'))
    up = up+1;
end

% appdir = strcat(ps(1:is(end-up)),'kernel/application/'); % application directory 
% libdir = strcat(ps(1:is(end-up)),'kernel/library/'); % library directory for linking
% valdir = strcat(ps(1:is(end-up)),'kernel/validation/main.cpp'); % validation directory for main cpp file
if (up==0)
    appdir = '/kernel/application/'; % application directory 
    libdir = '/library/'; % library directory for linking
    valdir = '/kernel/validation/main.cpp'; % validation directory for main cpp file        
elseif (up==1)
    appdir = '../kernel/application/'; % application directory 
    libdir = '../library/'; % library directory for linking
    valdir = '../kernel/validation/main.cpp'; % validation directory for main cpp file        
elseif (up==2)
    appdir = '../../kernel/application/'; % application directory 
    libdir = '../../library/'; % library directory for linking
    valdir = '../../kernel/validation/main.cpp'; % validation directory for main cpp file        
elseif (up==3)    
    appdir = '../../../kernel/application/'; % application directory 
    libdir = '../../../library/'; % library directory for linking
    valdir = '../../../kernel/validation/main.cpp'; % validation directory for main cpp file 
elseif (up==4)    
    appdir = '../../../../kernel/application/'; % application directory 
    libdir = '../../../../library/'; % library directory for linking
    valdir = '../../../../kernel/validation/main.cpp'; % validation directory for main cpp file        
elseif (up==5)    
    appdir = '../../../../../kernel/application/'; % application directory 
    libdir = '../../../../../library/'; % library directory for linking
    valdir = '../../../../../kernel/validation/main.cpp'; % validation directory for main cpp file            
elseif (up==6)    
    appdir = '../../../../../../kernel/application/'; % application directory 
    libdir = '../../../../../../library/'; % library directory for linking
    valdir = '../../../../../../kernel/validation/main.cpp'; % validation directory for main cpp file                
elseif (up==7)    
    appdir = '../../../../../../../kernel/application/'; % application directory 
    libdir = '../../../../../../../library/'; % library directory for linking
    valdir = '../../../../../../../kernel/validation/main.cpp'; % validation directory for main cpp file                    
elseif (up==8)    
    appdir = '../../../../../../../../kernel/application/'; % application directory 
    libdir = '../../../../../../../../library/'; % library directory for linking
    valdir = '../../../../../../../../kernel/validation/main.cpp'; % validation directory for main cpp file                        
end
appdir = strrep(appdir,'kernel',[version '/kernel']);
%libdir = strrep(libdir,'kernel',[version '/kernel']);
valdir = strrep(valdir,'kernel',[version '/kernel']);

fid  = fopen([appdir 'cpuTemplateApp.cpp'],'r');
f=fread(fid,'*char')';
fclose(fid);
f = strrep(f, 'ApplicationName', appname);    
fid  = fopen('cpuApp.cpp','w');
fprintf(fid,'%s',f);
fclose(fid);

fid  = fopen([appdir 'opuTemplateApp.cpp'],'r');
f=fread(fid,'*char')';
fclose(fid);
f = strrep(f, 'ApplicationName', appname);    
fid  = fopen('opuApp.cpp','w');
fprintf(fid,'%s',f);
fclose(fid);

fid  = fopen([appdir 'gpuTemplateApp.cu'],'r');
f=fread(fid,'*char')';
fclose(fid);
f = strrep(f, 'ApplicationName', appname);    
fid  = fopen('gpuApp.cu','w');
fprintf(fid,'%s',f);
fclose(fid);

if ismac
    if strcmp(blaslapack,'mkl')
        str = ['-L' libdir];
        blaslapack = strcat(str ,' -lmkl_intel_lp64 -lmkl_sequential -lmkl_core ');
        ompblaslapack = strcat(str ,' -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core ');
    end
        
    if ~isempty(cpucompiler)
    str = ['!' cpucompiler ' -fPIC -O3 -c cpuApp.cpp -fopenmp'];
    compilerstr{1} = str;
    eval(char(str));
    compilerstr{2} = str;
    str = '!ar -rvs cpuAppMac.a cpuApp.o';
    eval(char(str));
%     str = ['!' cpucompiler ' --shared -fPIC -O3 -fopenmp cpuApp.o -o libcpuAppMac.so'];
%     eval(char(str));

    str = ['!' cpucompiler ' -fPIC -O3 -c opuApp.cpp'];
    eval(char(str));
    compilerstr{3} = str;
    str = '!ar -rvs opuAppMac.a opuApp.o';
    compilerstr{4} = str;
    eval(char(str));
%     str = ['!' cpucompiler ' --shared -fPIC -O3 opuApp.o -o libopuAppMac.so'];
%     eval(char(str));        

    % compile serial code
    if ~isempty(blaslapack)
    disp('compile CPU C++ code');
    str = ['!' cpucompiler ' -std=c++11 ' cpuflags ' ' valdir ' -o cpu' appname ' '...
            libdir 'commonCoreMac.a ' libdir 'opuCoreMac.a opuAppMac.a -O2 -pthread -lm -ldl -Wl,-rpath,' libdir ' ' blaslapack];        
    compilerstr{5} = str;    
    eval(char(str));    
    end

    % compile OpenMP code
    if ~isempty(ompblaslapack)
    disp('compile OpenMP C++ code');
    str = ['!' cpucompiler ' -std=c++11 ' cpuflags ' ' valdir ' -o omp' appname ' '...
            libdir 'commonCoreMac.a ' libdir 'cpuCoreMac.a cpuAppMac.a -O2 -fopenmp -pthread -lm -ldl -Wl,-rpath,' libdir ' ' ompblaslapack];    
    eval(char(str));    
    compilerstr{6} = str;
    end    
    end        
    
    if ~isempty(mpicompiler)
    % compile MPI code
    if ~isempty(blaslapack)
    disp('compile MPI C++ code');     
    str = ['!' mpicompiler ' -std=c++11 -D _MPI ' cpuflags ' ' valdir ' -o mpi' appname ' '...
            libdir 'commonCoreMac.a ' libdir 'opuCoreMac.a opuAppMac.a -O2 -pthread -lm -ldl -Wl,-rpath,' libdir ' ' blaslapack];    
    eval(char(str));    
    compilerstr{7} = str;
    end
    
    % compile OpenMP-MPI code
    if ~isempty(ompblaslapack)
    disp('compile Open-MPI C++ code');         
    str = ['!' mpicompiler ' -std=c++11 -D _MPI ' cpuflags ' ' valdir ' -o ompmpi' appname ' '...
            libdir 'commonCoreMac.a ' libdir 'cpuCoreMac.a cpuAppMac.a -O2 -fopenmp -pthread -lm -ldl -Wl,-rpath,' libdir ' ' ompblaslapack];    
    eval(char(str));       
    compilerstr{8} = str;
    end 
    end
                
    if ~isempty(gpucompiler)
    str = ['!' gpucompiler ' -D_FORCE_INLINES -O3 -c --compiler-options ' char(39) '-fPIC' char(39) ' gpuApp.cu'];
    eval(char(str));
    compilerstr{9} = str;
    str = '!ar -rvs gpuAppMac.a gpuApp.o';
    eval(char(str));
    compilerstr{10} = str;
%     str = ['!' cpucompiler ' --shared -fPIC -O3 gpuApp.o -o libgpuAppMac.so'];
%     eval(char(str));    

    % compile GPU code
    if ~isempty(blaslapack) && ~isempty(cpucompiler)
    disp('compile GPU C++ code');         
    str = ['!' cpucompiler ' -std=c++11 -D _CUDA ' cpuflags ' ' valdir ' -o gpu' appname ' '...
            libdir 'commonCoreMac.a ' libdir 'opuCoreMac.a ' libdir 'gpuCoreMac.a' ' opuAppMac.a gpuAppMac.a '...
            '-O2 -pthread -lm -ldl ' gpuflags ' -Wl,-rpath,' libdir ' ' blaslapack];    
    eval(char(str));    
    compilerstr{11} = str;
    end
    
    % compile OpenMP-GPU code
    if ~isempty(ompblaslapack) && ~isempty(cpucompiler)
    disp('compile GPU-OpenMP C++ code');             
    str = ['!' cpucompiler ' -std=c++11 -D _CUDA ' cpuflags ' ' valdir ' -o gpuomp' appname ' '...
            libdir 'commonCoreMac.a ' libdir 'cpuCoreMac.a ' libdir 'gpuCoreMac.a' ' cpuAppMac.a gpuAppMac.a '...
            '-O2 -fopenmp -pthread -lm -ldl ' gpuflags ' -Wl,-rpath,' libdir ' ' ompblaslapack];    
    eval(char(str));  
    compilerstr{12} = str;
    end
    
    % compile GPU-MPI code
    if ~isempty(mpicompiler) && ~isempty(blaslapack)   
        disp('compile GPU-MPI C++ code');         
        str = ['!' mpicompiler ' -std=c++11 -D _MPI -D _CUDA ' cpuflags ' ' valdir ' -o gpumpi' appname ' '...
                libdir 'commonCoreMac.a ' libdir 'opuCoreMac.a ' libdir 'gpuCoreMac.a' ' opuAppMac.a gpuAppMac.a '...
                '-O2 -pthread -lm -ldl ' gpuflags ' -Wl,-rpath,' libdir ' ' blaslapack];    
        eval(char(str));            
        compilerstr{13} = str;
    end
    
    % compile GPU-OMP-MPI code
    if ~isempty(mpicompiler) && ~isempty(ompblaslapack)    
        disp('compile GPU-OpenMP-MPI C++ code');   
        str = ['!' mpicompiler ' -std=c++11 -D _MPI -D _CUDA ' cpuflags ' ' valdir ' -o gpuompmpi' appname ' '...
                libdir 'commonCoreMac.a ' libdir 'cpuCoreMac.a ' libdir 'gpuCoreMac.a' ' cpuAppMac.a gpuAppMac.a '...
                '-O2 -fopenmp -pthread -lm -ldl ' gpuflags ' -Wl,-rpath,' libdir ' ' ompblaslapack];    
        eval(char(str));            
        compilerstr{14} = str;
    end
    end    
elseif isunix
    if strcmp(blaslapack,'mkl')
        str = ['-L' libdir];
        blaslapack = strcat(str ,' -lmkl_intel_lp64 -lmkl_sequential -lmkl_core ');
        ompblaslapack = strcat(str ,' -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core ');
    end
    
    if ~isempty(cpucompiler)
    str = ['!' cpucompiler ' -fPIC -O3 -c cpuApp.cpp -fopenmp'];
    eval(char(str));
    compilerstr{1} = str;
    str = '!ar -rvs cpuAppLinux.a cpuApp.o';
    eval(char(str));
    compilerstr{2} = str;
%     str = ['!' cpucompiler ' --shared -fPIC -O3 -fopenmp cpuApp.o -o libcpuAppLinux.so'];
%     eval(char(str));

    str = ['!' cpucompiler ' -fPIC -O3 -c opuApp.cpp'];
    eval(char(str));
    compilerstr{3} = str;
    str = '!ar -rvs opuAppLinux.a opuApp.o';
    eval(char(str));
    compilerstr{4} = str;
%     str = ['!' cpucompiler ' --shared -fPIC -O3 opuApp.o -o libopuAppLinux.so'];
%     eval(char(str));        

    % compile serial code
    if ~isempty(blaslapack)
    disp('compile CPU C++ code');       
    str = ['!' cpucompiler ' -std=c++11 ' cpuflags ' ' valdir ' -o cpu' appname ' '... 
           libdir 'commonCoreLinux.a ' libdir 'opuCoreLinux.a opuAppLinux.a -O2 -pthread -lm -ldl -Wl,-rpath,'...
           libdir ' ' blaslapack];    
    eval(char(str));    
    compilerstr{5} = str;
    end

    % compile OpenMP code
    if ~isempty(ompblaslapack)
    disp('compile OpenMP C++ code');           
    str = ['!' cpucompiler ' -std=c++11 ' cpuflags ' ' valdir ' -o omp' appname ' '...
           libdir 'commonCoreLinux.a ' libdir 'cpuCoreLinux.a cpuAppLinux.a -O2 -fopenmp -pthread -lm -ldl -Wl,-rpath,'...
           libdir ' ' ompblaslapack];    
    eval(char(str));        
    compilerstr{6} = str;
    end    
    end        
    
    if ~isempty(mpicompiler)
    % compile MPI code
    if ~isempty(blaslapack)
    disp('compile MPI C++ code');        
    str = ['!' mpicompiler ' -std=c++11 -D _MPI ' cpuflags ' ' valdir ' -o mpi' appname ' '...
            libdir 'commonCoreLinux.a ' libdir 'opuCoreLinux.a opuAppLinux.a -O2 -pthread -lm -ldl -Wl,-rpath,'...
            libdir ' ' blaslapack];    
    eval(char(str));       
    compilerstr{7} = str;
    end
   
    % compile OpenMP-MPI code
    if ~isempty(ompblaslapack)
    disp('compile OpenMP-MPI C++ code');            
    str = ['!' mpicompiler ' -std=c++11 -D _MPI ' cpuflags ' ' valdir ' -o ompmpi' appname ' '...
            libdir 'commonCoreLinux.a ' libdir 'cpuCoreLinux.a cpuAppLinux.a -O2 -fopenmp -pthread -lm -ldl -Wl,-rpath,'...
            libdir ' ' ompblaslapack];    
    eval(char(str));    
    compilerstr{8} = str;
    end
    end
                
    if ~isempty(gpucompiler)
    str = ['!' gpucompiler ' -D_FORCE_INLINES -O3 -c --compiler-options ' char(39) '-fPIC' char(39) ' gpuApp.cu'];
    eval(char(str));
    compilerstr{9} = str;
    str = '!ar -rvs gpuAppLinux.a gpuApp.o';
    eval(char(str));
    compilerstr{10} = str;
%     str = ['!' cpucompiler ' --shared -fPIC -O3 gpuApp.o -o libgpuAppLinux.so'];
%     eval(char(str));    

    % compile GPU code
    if ~isempty(blaslapack) && ~isempty(cpucompiler)
    disp('compile GPU C++ code');             
    str = ['!' cpucompiler ' -std=c++11 -D _CUDA ' cpuflags ' ' valdir ' -o gpu' appname ' '...
            libdir 'commonCoreLinux.a ' libdir 'opuCoreLinux.a ' libdir 'gpuCoreLinux.a' ' opuAppLinux.a gpuAppLinux.a '...
            '-O2 -pthread -lm -ldl ' gpuflags ' -Wl,-rpath,' libdir ' ' blaslapack];    
    eval(char(str));    
    compilerstr{11} = str;
    end
    
    % compile OpenMP-GPU code
    if ~isempty(ompblaslapack) && ~isempty(cpucompiler)
    disp('compile GPU-OpenMP C++ code');                 
    str = ['!' cpucompiler ' -std=c++11 -D _CUDA ' cpuflags ' ' valdir ' -o gpuomp' appname ' '...
            libdir 'commonCoreLinux.a ' libdir 'cpuCoreLinux.a ' libdir 'gpuCoreLinux.a' ' cpuAppLinux.a gpuAppLinux.a '...
            '-O2 -fopenmp -pthread -lm -ldl ' gpuflags ' -Wl,-rpath,' libdir ' ' ompblaslapack];    
    eval(char(str));  
    compilerstr{12} = str;
    end
    
    % compile GPU-MPI code
    if ~isempty(mpicompiler) && ~isempty(blaslapack)        
        disp('compile GPU-MPI C++ code');                 
        str = ['!' mpicompiler ' -std=c++11 -D _MPI -D _CUDA ' cpuflags ' ' valdir ' -o gpumpi' appname ' '...
                libdir 'commonCoreLinux.a ' libdir 'opuCoreLinux.a ' libdir 'gpuCoreLinux.a' ' opuAppLinux.a gpuAppLinux.a '...
                '-O2 -pthread -lm -ldl ' gpuflags ' -Wl,-rpath,' libdir ' ' blaslapack];    
        eval(char(str));      
        compilerstr{13} = str;
    end
    
    % compile GPU-OMP-MPI code
    if ~isempty(mpicompiler) && ~isempty(ompblaslapack)        
        disp('compile GPU-OpenMP-MPI C++ code');                 
        str = ['!' mpicompiler ' -std=c++11 -D _MPI -D _CUDA ' cpuflags ' ' valdir ' -o gpuompmpi' appname ' '...
                libdir 'commonCoreLinux.a ' libdir 'cpuCoreLinux.a ' libdir 'gpuCoreLinux.a' ' cpuAppLinux.a gpuAppLinux.a '...
                '-O2 -fopenmp -pthread -lm -ldl ' gpuflags ' -Wl,-rpath,' libdir ' ' ompblaslapack];    
        eval(char(str));          
        compilerstr{14} = str;
    end
    end
end

if deletesourcefiles==1
    delete *.o;
    delete *.a;    
    delete *.cpp;
    delete *.cu;
end


