import shutil
import os
from sys import platform
from genlib import genlib

def checkcompilers(app):

    cpustatus = shutil.which(app['cpucompiler']);
    if cpustatus != None:
        print("Using " + app['cpucompiler'] + " compiler for CPU source code");
    else:
        print(app['cpucompiler'] + " compiler is not available on your system.");
        cpustatus = shutil.which("g++");
        if cpustatus != None:
            print("However, g++ compiler is available on your system. It will be used for compiling CPU source code.");
            app['cpucompiler'] = "g++";
        else:
            error("C++ compiler is not available on your system.");

    if app['mpiprocs']>1:
        mpistatus = shutil.which(app['mpicompiler']);
        if mpistatus != None:
            print("Using " + app['mpicompiler'] + " compiler for MPI source code");
        else:
            print(app['mpicompiler'] + " compiler is not available on your system.");
            mpistatus = shutil.which("mpicxx");
            if mpistatus != None:
                print("However, mpicxx compiler is available on your system. It will be used for compiling MPI source code.");
                app['mpicompiler'] = "mpicxx";
                app['mpirun'] = "mpirun";
            else:
                error("MPI compiler is not available on your system.");

    if app['platform'] == "gpu":
        gpustatus = shutil.which(app['gpucompiler']);
        if gpustatus != None:
            print("Using " + app['gpucompiler'] + " compiler for GPU source code");
        else:
            print(app['gpucompiler'] + " compiler is not available on your system.");
            gpustatus = shutil.which("nvcc");
            if gpustatus != None:
                print("However, nvcc compiler is available on your system. It will be used for compiling GPU source code.");
                app['gpucompiler'] = "nvcc";
            else:
                error("GPU compiler is not available on your system.");

    mdir = os.getcwd();
    ii = mdir.find("Exasim");
    coredir = mdir[0:(ii+6)] + "/lib";

    cpulib = 0;
    gpulib = 0;
    if platform == "darwin":
        if (os.path.isfile(coredir + "/Mac/commonCore.a")) and (os.path.isfile(coredir + "/Mac/opuCore.a")):
            cpulib = 1;
        if (os.path.isfile(coredir + "/Mac/gpuCore.a")):
            gpulib = 1;
    elif platform == "linux" or platform == "linux2":
        if (os.path.isfile(coredir + "/Linux/commonCore.a")) and (os.path.isfile(coredir + "/Linux/opuCore.a")):
            cpulib = 1;
        if (os.path.isfile(coredir + "/Linux/gpuCore.a")):
            gpulib = 1;
    elif platform == "win32":
        if (os.path.isfile(coredir + "/Windows/commonCore.a")) and (os.path.isfile(coredir + "/Windows/opuCore.a")):
            cpulib = 1;
        if (os.path.isfile(coredir + "/Windows/gpuCore.a")):
            gpulib = 1;

    if cpulib==0:
        print("Generating CPU core libraries.");
        genlib(app['cpucompiler'], "", coredir);

    if gpulib==0:
        if (len(app['gpucompiler'])>0) and (app['platform'] == "gpu"):
            print("Generating GPU core library.");
            genlib("", app['gpucompiler'], coredir);

    return cpulib, gpulib
