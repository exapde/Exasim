from sys import platform
import shutil, os, sys
import subprocess

version = "src";
cdir = os.getcwd(); ii = cdir.find("Exasim");
versiondir = cdir[0:(ii+6)] + "/"  + version + "/Python";

#sys.path.append(cdir[0:(ii+6)] + "/Installation");
from findinstallexec import findinstallexec

if platform == "darwin":
    brewstatus0 = shutil.which("brew");
    brewstatus1 = shutil.which("/usr/bin/brew");
    brewstatus2 = shutil.which("/usr/local/bin/brew");
    if brewstatus0 != None:
        brew = "brew";
    elif brewstatus1 != None:
        brew = "/usr/bin/brew";
    elif brewstatus2 != None:
        brew = "/usr/local/bin/brew";
    else:
        print("Homebrew is not installed on your computer.\n")
        print("Please visit https://brew.sh/ to install Homebrew.\n");
        print("Run this script again after installing Homebrew.\n");
        print("Installation terminated because Homebrew is not installed yet.\n");
        raise ValueError()
else:
    brew = "";

print("Exasim find and install the required packages...\n");

gcc, appinstall = findinstallexec("g++", "gcc", brew, 0);
if appinstall==1:
    gcc = findinstallexec("g++", "gcc", brew, 0);

mpi, appinstall = findinstallexec("mpicxx", "openmpi", brew, 0);
if appinstall==1:
    mpi = findinstallexec("mpicxx", "openmpi", brew, 0);

nvcc, appinstall = findinstallexec("nvcc", "nvcc", brew, 10);
if appinstall==1:
    nvcc = findinstallexec("nvcc", "nvcc", brew, 10);

metis, appinstall = findinstallexec("mpmetis", "metis", brew, 0);
if appinstall==1:
    metis = findinstallexec("mpmetis", "metis", brew, 0);

#python = findinstallexec("python", "python", brew, 0);
gmsh, appinstall = findinstallexec("gmsh", "gmsh", brew, 0);
if appinstall==1:
    gmsh = findinstallexec("gmsh", "gmsh", brew, 0);

if platform == "darwin":
    #julia = findinstallexec("julia", "julia", brew, 1);
    paraview,appinstall = findinstallexec("paraview", "paraview", brew, 1);
    if appinstall==1:
        paraview = findinstallexec("paraview", "paraview", brew, 1);
else:
    #julia = findinstallexec("julia", "julia", brew, 0);
    paraview,appinstall = findinstallexec("paraview", "paraview", brew, 0);
    if appinstall==1:
        paraview = findinstallexec("paraview", "paraview", brew, 0);

if platform == "darwin":
    print("Installing " +  "openblas" +  " via brew.");
    os.system(brew +  " install openblas");
    print("Installing " +  "lapack" +  " via brew.");
    os.system(brew +  " install lapack");
else:
    print("Installing " +  "blas" +  " via apt.");
    os.system("sudo apt install libblas-dev");
    print("Installing " +  "lapack" +  " via apt.");
    os.system("sudo apt install liblapack-dev");

filem = versiondir +  "/Preprocessing/initializepde.py";
fid = open(filem, "r");
text = fid.read();
fid.close();

q = """ " """;
q = q[1];

oldgcc = "pde['cpucompiler'] = " +  q +  "g++" +  q;
p = q +  gcc +  q;
newgcc = "pde['cpucompiler']  = " +  p;
newtext = text.replace(oldgcc, newgcc);

oldmpi = "pde['mpicompiler'] = " +  q +  "mpicxx" +  q;
p = q +  mpi +  q;
newmpi = "pde['mpicompiler'] = " +  p;
newtext = newtext.replace(oldmpi, newmpi);

oldnvcc = "pde['gpucompiler'] = " +  q +  "nvcc" +  q;
p = q +  nvcc +  q;
newnvcc = "pde['gpucompiler'] = " +  p;
newtext = newtext.replace(oldnvcc, newnvcc);

oldmpirun = "pde['mpirun'] = " +  q +  "mpirun" +  q;
p = q +  mpi +  q;
p = p.replace("mpicxx", "mpirun");
newmpirun = "pde['mpirun'] = " +  p;
newtext = newtext.replace(oldmpirun, newmpirun);

oldmetis = "pde['metis'] = " +  q +  "mpmetis" +  q;
p = q +  metis +  q;
newmetis = "pde['metis'] = " +  p;
newtext = newtext.replace(oldmetis, newmetis);

oldgmsh = "pde['gmsh'] = " +  q +  "gmsh" +  q;
p = q +  gmsh +  q;
newgmsh = "pde['gmsh'] = " +  p;
newtext = newtext.replace(oldgmsh, newgmsh);

oldparaview = "pde['paraview'] = " +  q +  "paraview" +  q;
p = q +  paraview +  q;
newparaview = "pde['paraview'] = " +  p;
newtext = newtext.replace(oldparaview, newparaview);

fid = open(filem, "w");
fid.write(newtext);
fid.close();
