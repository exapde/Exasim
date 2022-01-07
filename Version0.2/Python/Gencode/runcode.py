import os

def runcode(app, numpde):

    print("run code...");

    if numpde==1:
        mystr = app['appname'] + " 1 datain/ dataout/out";
    else:
        mystr = app['appname'] + " " + str(numpde);
        for i in range(1, numpde+1):
            mystr = mystr + " datain" + str(i) + "/ " + "dataout" + str(i) + "/out";

    mpirun = app['mpirun'];
    if app['platform'] == "cpu":
        if app['mpiprocs']==1:
            mystr = "./app/serial" + mystr;                    
        else:
            mystr = mpirun + " -np " + str(app['mpiprocs']) + " ./app/mpi" + mystr;        
        os.system(mystr);
    elif app['platform'] == "gpu":
        if app['mpiprocs']==1:
            mystr = "./app/gpu" + mystr;   
        else:
            mystr = mpirun + " -np " + str(app['mpiprocs']) + " ./app/gpumpi" + mystr;        
        os.system(mystr);
    else:
        error("app['platform'] must be either cpu or gpu");

    return mystr
