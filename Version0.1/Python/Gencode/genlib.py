import os
from sys import platform

def genlib(cpucompiler, gpucompiler, coredir):

    mydir = os.getcwd();
    os.chdir(coredir);

    if len(cpucompiler)>0:
        str = cpucompiler + " -fPIC -O3 -c commonCore.cpp";
        os.system(str);
        str = "ar rvs commonCore.a commonCore.o";
        os.system(str);

        str = cpucompiler + " -fPIC -O3 -c opuCore.cpp";
        os.system(str);
        str = "ar rvs opuCore.a opuCore.o";
        os.system(str);

    if len(gpucompiler)>0:
        str = gpucompiler + " -D_FORCE_INLINES -O3 -c --compiler-options '-fPIC' gpuCore.cu";
        os.system(str);
        str = "ar rvs gpuCore.a gpuCore.o";
        os.system(str);

    if platform == "darwin":
    #    mv("*.a",pwd() * "/Mac");
        if os.path.isdir("Mac")==False:
            os.mkdir("Mac")

        if os.path.isfile("commonCore.a"):
            os.system("mv commonCore.a Mac");

        if os.path.isfile("opuCore.a"):
            os.system("mv opuCore.a Mac");

        if os.path.isfile("cpuCore.a"):
            os.system("mv cpuCore.a Mac");

        if os.path.isfile("gpuCore.a"):
            os.system("mv gpuCore.a Mac");
    elif platform == "linux" or platform == "linux2":
        if os.path.isdir("Linux")==False:
            os.mkdir("Linux")

        if os.path.isfile("commonCore.a"):
            os.system("mv commonCore.a Linux");

        if os.path.isfile("opuCore.a"):
            os.system("mv opuCore.a Linux");

        if os.path.isfile("cpuCore.a"):
            os.system("mv cpuCore.a Linux");

        if os.path.isfile("gpuCore.a"):
            os.system("mv gpuCore.a Linux");

    elif platform == "win32":
        if os.path.isdir("Windows")==False:
            os.mkdir("Windows")

        if os.path.isfile("commonCore.a"):
            os.system("mv commonCore.a Windows");

        if os.path.isfile("opuCore.a"):
            os.system("mv opuCore.a Windows");

        if os.path.isfile("cpuCore.a"):
            os.system("mv cpuCore.a Windows");

        if os.path.isfile("gpuCore.a"):
            os.system("mv gpuCore.a Windows");

    # os.remove("commonCore.a");
    # os.remove("opuCore.a");
    # os.remove("cpuCore.a");
    # os.remove("gpuCore.a");
    # os.remove("commonCore.o");
    # os.remove("opuCore.o");
    # os.remove("cpuCore.o");
    # os.remove("gpuCore.o");

    os.chdir(mydir);

    return 0
