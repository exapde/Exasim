import os

def runcode(app):

    print("run code...");
    mpirun = app['mpirun'];
    if app['platform'] == "cpu":
        if app['mpiprocs']==1:
            mystr = "./app/serial" + app['appname'] + " datain/ dataout/out";
        else:
            mystr = mpirun + " -np " + str(app['mpiprocs']) + " ./app/mpi" + app['appname'] + " datain/ dataout/out";

        os.system(mystr);
    elif app['platform'] == "gpu":
        if app['mpiprocs']==1:
            mystr = "./app/gpu" + app['appname'] + " datain/ dataout/out";
        else:
            mystr = mpirun + " -np " + str(app['mpiprocs']) + " ./app/gpumpi" + app['appname'] + " datain/ dataout/out";
        os.system(mystr);
    else:
        error("app['platform'] must be either cpu or gpu");

    return mystr
