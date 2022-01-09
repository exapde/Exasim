import shutil
import os
from sys import platform, exit
from genlib import genlib

def setcompilers(app):

    cpustatus0 = shutil.which(app['cpucompiler']);
    cpustatus1 = shutil.which("g++");
    cpustatus2 = shutil.which("/usr/bin/g++");
    cpustatus3 = shutil.which("/usr/local/bin/g++");
    cpustatus4 = shutil.which("/opt/local/bin/g++");

    if cpustatus0 != None:
        print("Using " + app['cpucompiler'] + " compiler for CPU source code");
    elif cpustatus1 != None:
        app['cpucompiler'] = "g++";
        print("Using g++ compiler for CPU source code");
    elif cpustatus2 != None:
        app['cpucompiler'] = "/usr/bin/g++";
        print("Using /usr/bin/g++ compiler for CPU source code");
    elif cpustatus3 != None:
        app['cpucompiler'] = "/usr/local/bin/g++";
        print("Using /usr/local/bin/g++ compiler for CPU source code");
    elif cpustatus4 != None:
        app['cpucompiler'] = "/opt/local/bin/g++";
        print("Using /opt/local/bin/g++ compiler for CPU source code");
    else:
        exit("Exasim search in /usr/bin, /usr/local/bin, and /opt/local/bin and could not find C++ compiler. Please see the documentation to install it. After installation, please set its path to app['cpucompiler']");

    if app['mpiprocs']>1:
        mpistatus0 = shutil.which(app['mpicompiler']);
        mpistatus1 = shutil.which("mpicxx");
        mpistatus2 = shutil.which("/usr/bin/mpicxx");
        mpistatus3 = shutil.which("/usr/local/bin/mpicxx");
        mpistatus4 = shutil.which("/opt/local/bin/mpicxx");
        mpistatus5 = shutil.which("mpicxx-openmpi-mp");
        mpistatus6 = shutil.which("/usr/bin/mpicxx-openmpi-mp");
        mpistatus7 = shutil.which("/usr/local/bin/mpicxx-openmpi-mp");
        mpistatus8 = shutil.which("/opt/local/bin/mpicxx-openmpi-mp");
        if mpistatus0 != None:
            print("Using " + app['mpicompiler'] + " compiler for MPI source code");
        elif mpistatus5 != None:
            app['mpicompiler'] = "mpicxx-openmpi-mp";
            app['mpirun'] = "mpirun-openmpi-mp";
            print("Using mpicxx-openmpi-mp compiler for MPI source code");
        elif mpistatus6 != None:
            app['mpicompiler'] = "/usr/bin/mpicxx-openmpi-mp";
            app['mpirun'] = "/usr/bin/mpirun-openmpi-mp";
            print("Using /usr/bin/mpicxx-openmpi-mp compiler for MPI source code");
        elif mpistatus7 != None:
            app['mpicompiler'] = "/usr/local/bin/mpicxx-openmpi-mp";
            app['mpirun'] = "/usr/local/bin/mpirun-openmpi-mp";
            print("Using /usr/local/bin/mpicxx-openmpi-mp compiler for MPI source code");
        elif mpistatus8 != None:
            app['mpicompiler'] = "/opt/local/bin/mpicxx-openmpi-mp";
            app['mpirun'] = "/opt/local/bin/mpirun-openmpi-mp";
            print("Using opt/local/bin/mpicxx-openmpi-mp compiler for MPI source code");
        elif mpistatus1 != None:
            app['mpicompiler'] = "mpicxx";
            app['mpirun'] = "mpirun";
            print("Using mpicxx compiler for MPI source code");
        elif mpistatus2 != None:
            app['mpicompiler'] = "/usr/bin/mpicxx";
            app['mpirun'] = "/usr/bin/mpirun";
            print("Using /usr/bin/mpicxx compiler for MPI source code");
        elif mpistatus3 != None:
            app['mpicompiler'] = "/usr/local/bin/mpicxx";
            app['mpirun'] = "/usr/local/bin/mpirun";
            print("Using /usr/local/bin/mpicxx compiler for MPI source code");
        elif mpistatus4 != None:
            app['mpicompiler'] = "/opt/local/bin/mpicxx";
            app['mpirun'] = "/opt/local/bin/mpirun";
            print("Using opt/local/bin/mpicxx compiler for MPI source code");
        else:
            exit("Exasim search in /usr/bin, /usr/local/bin, and /opt/local/bin and could not find MPI compiler. Please see the documentation to install it. After installation, please set its path to app['mpicompiler']");
        
        app['mpirun'] = app['mpicompiler'].replace("mpicxx", "mpirun");

    if app['platform'] == "gpu":
        gpustatus0 = shutil.which(app['gpucompiler']);
        gpustatus1 = shutil.which("nvcc");
        gpustatus2 = shutil.which("/usr/bin/nvcc");
        gpustatus3 = shutil.which("/usr/local/bin/nvcc");
        gpustatus4 = shutil.which("/opt/local/bin/nvcc");

        if gpustatus0 != None:
            print("Using " + app['gpucompiler'] + " compiler for GPU source code");
        elif gpustatus1 != None:
            app['gpucompiler'] = "nvcc";
            print("Using nvcc compiler for GPU source code");
        elif gpustatus2 != None:
            app['gpucompiler'] = "/usr/bin/nvcc";
            print("Using /usr/bin/nvcc compiler for GPU source code");
        elif gpustatus3 != None:
            app['gpucompiler'] = "/usr/local/bin/nvcc";
            print("Using /usr/local/bin/nvcc compiler for GPU source code");
        elif gpustatus4 != None:
            app['gpucompiler'] = "opt/local/bin/nvcc";
            print("Using opt/local/bin/nvcc compiler for GPU source code");
        else:
            exit("Exasim search in /usr/bin, /usr/local/bin, and /opt/local/bin and could not find NVCC compiler. Please see the documentation to install it. After installation, please set its path to app['gpucompiler']");

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

    return app
