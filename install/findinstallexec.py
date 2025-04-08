from sys import platform
import shutil, os
import subprocess

def findinstallexec(filename, appname, brew, searchopt):

    print("Finding " +  appname +  "...\n");

    dirname = "";
    appinstall = 0;

    status = shutil.which(filename);
    if status != None:
        dirname = filename;
        print("Exasim found " +  dirname +  "\n");
        return dirname, appinstall;

    status = shutil.which("/usr/bin/" + filename);
    if status != None:
        dirname = "/usr/bin/" +  filename;
        print("Exasim found " +  dirname +  "\n");
        return dirname, appinstall;

    status = shutil.which("/usr/local/bin/" + filename);
    if status != None:
        dirname = "/usr/local/bin/" +  filename;
        print("Exasim found " +  dirname +  "\n");
        return dirname, appinstall;

    status = shutil.which("/opt/local/bin/" + filename);
    if status != None:
        dirname = "/opt/local/bin/" +  filename;
        print("Exasim found " +  dirname +  "\n");
        return dirname, appinstall;

    if searchopt==1:
        if platform == "darwin":
            result = subprocess.run(['find', '/Applications', '-name', filename, '-type', 'f'], stdout=subprocess.PIPE);
            a = result.stdout.decode('utf-8');
            if len(a)>0:
                ii = a.find("/MacOS/" + filename);
                dirname = a[0:ii] + "/MacOS/" + filename;
                print("Exasim found " +  dirname +  "\n");
                return dirname, appinstall;

    print("Exasim could not find " +  appname +  " on your computer.\n");
    if searchopt==10:
        dirname = filename;
        print("CUDA Toolkit is not found on your system.\n");
        print("If you have Nividia GPUs on your system, please visit https://docs.nvidia.com/cuda/ to install it.\n");
        return dirname, appinstall;
    else:
        appinstall = 1;
        if platform == "darwin":
            print("Installing " +  appname +  " via brew.\n");
            if searchopt==1:
                os.systemn(brew +  "install --cask " +  appname);
            else:
                os.system(brew +  " install " +  appname);

        elif platform == "linux" or platform == "linux2":
            print("Installing " +  appname +  " via apt.\n");
            os.system("sudo apt install " +  appname);
        elif platform == "win32":
            print("Installing " +  appname +  " via apt.\n");
            os.system("sudo apt install " +  appname);

        return dirname, appinstall;
