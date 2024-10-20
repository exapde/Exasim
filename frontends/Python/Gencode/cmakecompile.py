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
            comstr = ["cmake", "-D", "EXASIM_NOMPI=ON", "-D", "EXASIM_MPI=OFF", "-D", "EXASIM_CUDA=ON", "../install"]
        else:
            comstr = ["cmake", "-D", "EXASIM_NOMPI=ON", "-D", "EXASIM_MPI=OFF", "-D", "EXASIM_CUDA=OFF", "../install"]
    else:
        if pde['platform'] == "gpu":
            comstr = ["cmake", "-D", "EXASIM_NOMPI=OFF", "-D", "EXASIM_MPI=ON", "-D", "EXASIM_CUDA=ON", "../install"]
        else:
            comstr = ["cmake", "-D", "EXASIM_NOMPI=OFF", "-D", "EXASIM_MPI=ON", "-D", "EXASIM_CUDA=OFF", "../install"]

    subprocess.run(comstr)
    subprocess.run(["cmake", "--build", "."])

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
