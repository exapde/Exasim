if ismac
    !/opt/local/bin/g++-mp-7 -fPIC -O3 -c commonCore.cpp
    !/opt/local/bin/g++-mp-7 -dynamiclib -fPIC -o libcommonCoreMac.dylib commonCore.o
    !ar rvs commonCoreMac.a commonCore.o

    !/opt/local/bin/g++-mp-7 -fPIC -O3 -c opuCore.cpp
    !/opt/local/bin/g++-mp-7 -dynamiclib -fPIC -o libopuCoreMac.dylib opuCore.o
    !ar rvs opuCoreMac.a opuCore.o

    !/opt/local/bin/g++-mp-7 -fPIC -O3 -c cpuCore.cpp -fopenmp
    !/opt/local/bin/g++-mp-7 -dynamiclib -fPIC -o libcpuCoreMac.dylib cpuCore.o -fopenmp
    !ar rvs cpuCoreMac.a cpuCore.o     
    
    copyfile('*.a','../../../library');
    copyfile('*.dylib','../../../library');
elseif isunix
    !g++ -fPIC -O3 -c commonCore.cpp
    !g++ --shared commonCore.o -o libcommonCoreLinux.so
    !ar rvs commonCoreLinux.a commonCore.o    

    !g++ -fPIC -O3 -c opuCore.cpp
    !g++ --shared opuCore.o -o libopuCoreLinux.so
    !ar rvs opuCoreLinux.a opuCore.o    

    !g++ -fPIC -O3 -c cpuCore.cpp -fopenmp
    !g++ --shared cpuCore.o -o libcpuCoreLinux.so
    !ar rvs cpuCoreLinux.a cpuCore.o        
    
    %!nvcc -D_FORCE_INLINES -O3 -c gpuCore.cu
    !nvcc -D_FORCE_INLINES -O3 -c --compiler-options '-fPIC' gpuCore.cu
    !g++ --shared gpuCore.o -o libgpuCoreLinux.so    
    !ar -rvs gpuCoreLinux.a gpuCore.o    
    
    copyfile('*.a','../../../library');
    copyfile('*.so','../../../library');
end

