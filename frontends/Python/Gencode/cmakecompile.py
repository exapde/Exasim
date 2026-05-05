import os
import subprocess

def cmakecompile(pde):
    
    print("Compile C++ Exasim code using cmake...")

    cdir = os.getcwd()    
    os.chdir(pde['buildpath'])

    files_to_delete = ["cpuEXASIM", "cpumpiEXASIM", "gpuEXASIM", "gpumpiEXASIM"]
    for file in files_to_delete:
        if os.path.exists(file):
            os.remove(file)

    if pde['mpiprocs'] == 1:
        if pde['platform'] == "gpu":
            comstr = ["cmake", "-D", "EXASIM_NOMPI=ON", "-D", "EXASIM_MPI=OFF", "-D", "EXASIM_CUDA=ON", "-D", "EXASIM_BUILD_LIBRARY_EXAMPLES=OFF", "../install"]
            target = "gpuEXASIM"
        elif pde['platform'] == "hip":
            comstr = ["cmake", "-D", "CMAKE_CXX_COMPILER=hipcc", "-D", "EXASIM_NOMPI=ON", "-D", "EXASIM_HIP=ON", "-D", "EXASIM_BUILD_LIBRARY_EXAMPLES=OFF", "../install"]
            target = "gpuEXASIM"
        else:
            comstr = ["cmake", "-D", "EXASIM_NOMPI=ON", "-D", "EXASIM_MPI=OFF", "-D", "EXASIM_CUDA=OFF", "-D", "EXASIM_BUILD_LIBRARY_EXAMPLES=OFF", "../install"]
            target = "cpuEXASIM"
    else:
        if pde['platform'] == "gpu":
            comstr = ["cmake", "-D", "EXASIM_NOMPI=OFF", "-D", "EXASIM_MPI=ON", "-D", "EXASIM_CUDA=ON", "-D", "EXASIM_BUILD_LIBRARY_EXAMPLES=OFF", "../install"]
            target = "gpumpiEXASIM"
        elif pde['platform'] == "hip":
            comstr = ["cmake", "-D", "CMAKE_CXX_COMPILER=hipcc", "-D", "EXASIM_NOMPI=OFF", "-D", "EXASIM_MPI=ON", "-D", "EXASIM_HIP=ON", "-D", "EXASIM_BUILD_LIBRARY_EXAMPLES=OFF", "../install"]
            target = "gpumpiEXASIM"
        else:
            comstr = ["cmake", "-D", "EXASIM_NOMPI=OFF", "-D", "EXASIM_MPI=ON", "-D", "EXASIM_CUDA=OFF", "-D", "EXASIM_BUILD_LIBRARY_EXAMPLES=OFF", "../install"]
            target = "cpumpiEXASIM"

    # check=True so configure / build failures surface here, not as
    # confusing downstream errors when the runner can't find the
    # binary or fetchsolution sees missing dataout files.
    subprocess.run(comstr, check=True)
    subprocess.run(["cmake", "--build", ".", "--target", target], check=True)

    os.chdir(cdir)

    return comstr

    # print("Run C++ Exasim code ...")
    # pdenum = " " + str(numpde) + " "
    # mpirun = pde['mpirun']
    # DataPath = pde['buildpath']    

    # if pde['platform'] == "cpu":
    #     if pde['mpiprocs'] == 1:
    #         cmd = f"./cpuEXASIM {pdenum} {DataPath}/datain/ {DataPath}/dataout/out"
    #     else:
    #         cmd = f"{mpirun} -np {pde['mpiprocs']} ./cpumpiEXASIM {pdenum} {DataPath}/datain/ {DataPath}/dataout/out"
    # elif pde['platform'] == "gpu":
    #     if pde['mpiprocs'] == 1:
    #         cmd = f"./gpuEXASIM {pdenum} {DataPath}/datain/ {DataPath}/dataout/out"
    #     else:
    #         cmd = f"{mpirun} -np {pde['mpiprocs']} ./gpumpiEXASIM {pdenum} {DataPath}/datain/ {DataPath}/dataout/out"

    # start_time = time.time()    
    # subprocess.run(cmd, shell=True)    
    # end_time = time.time()
    # print(f"Elapsed time: {end_time - start_time} seconds")

    # os.chdir(cdir)
    
    # return cmd
