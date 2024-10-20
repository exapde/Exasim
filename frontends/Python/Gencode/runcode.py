import os
import subprocess
import time

def runcode(pde, numpde):

    print("Run C++ Exasim code ...")

    cdir = os.getcwd()    
    os.chdir(pde['buildpath'])

    pdenum = " " + str(numpde) + " "
    mpirun = pde['mpirun']
    DataPath = pde['buildpath']

    # buildpath = pde['buildpath']    
    # if pde['platform'] == "cpu":
    #     if pde['mpiprocs'] == 1:
    #         cmd = f"{buildpath}/cpuEXASIM {pdenum} {DataPath}/datain/ {DataPath}/dataout/out"
    #     else:
    #         cmd = f"{mpirun} -np {pde['mpiprocs']} {buildpath}/cpumpiEXASIM {pdenum} {DataPath}/datain/ {DataPath}/dataout/out"
    # elif pde['platform'] == "gpu":
    #     if pde['mpiprocs'] == 1:
    #         cmd = f"{buildpath}/gpuEXASIM {pdenum} {DataPath}/datain/ {DataPath}/dataout/out"
    #     else:
    #         cmd = f"{mpirun} -np {pde['mpiprocs']} {buildpath}/gpumpiEXASIM {pdenum} {DataPath}/datain/ {DataPath}/dataout/out"

    if pde['platform'] == "cpu":
        if pde['mpiprocs'] == 1:
            cmd = f"./cpuEXASIM {pdenum} {DataPath}/datain/ {DataPath}/dataout/out"
        else:
            cmd = f"{mpirun} -np {pde['mpiprocs']} ./cpumpiEXASIM {pdenum} {DataPath}/datain/ {DataPath}/dataout/out"
    elif pde['platform'] == "gpu":
        if pde['mpiprocs'] == 1:
            cmd = f"./gpuEXASIM {pdenum} {DataPath}/datain/ {DataPath}/dataout/out"
        else:
            cmd = f"{mpirun} -np {pde['mpiprocs']} ./gpumpiEXASIM {pdenum} {DataPath}/datain/ {DataPath}/dataout/out"

    start_time = time.time()    
    subprocess.run(cmd, shell=True)    
    end_time = time.time()
    
    print(f"Elapsed time: {end_time - start_time} seconds")
    
    os.chdir(cdir)

    return cmd

    # print("run code...");

    # if numpde==1:
    #     mystr = app['appname'] + " 1 datain/ dataout/out";
    # else:
    #     mystr = app['appname'] + " " + str(numpde);
    #     for i in range(1, numpde+1):
    #         mystr = mystr + " datain" + str(i) + "/ " + "dataout" + str(i) + "/out";

    # mpirun = app['mpirun'];
    # if app['platform'] == "cpu":
    #     if app['mpiprocs']==1:
    #         mystr = "./app/serial" + mystr;                    
    #     else:
    #         mystr = mpirun + " -np " + str(app['mpiprocs']) + " ./app/mpi" + mystr;        
    #     os.system(mystr);
    # elif app['platform'] == "gpu":
    #     if app['mpiprocs']==1:
    #         mystr = "./app/gpu" + mystr;   
    #     else:
    #         mystr = mpirun + " -np " + str(app['mpiprocs']) + " ./app/gpumpi" + mystr;        
    #     os.system(mystr);
    # else:
    #     error("app['platform'] must be either cpu or gpu");

    # return mystr
