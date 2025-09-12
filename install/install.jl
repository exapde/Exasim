include("findinstallexec.jl");

version = "src";

if Sys.isapple()
    brewstatus0 = Sys.which("brew");
    brewstatus1 = Sys.which("/usr/bin/brew");
    brewstatus2 = Sys.which("/usr/local/bin/brew");
    if brewstatus0 != nothing
        brew = "brew";
    elseif brewstatus1 != nothing
        brew = "/usr/bin/brew";
    elseif brewstatus2 != nothing
        brew = "/usr/local/bin/brew";
    else
        print("Homebrew is not installed on your computer.\n")
        print("Please visit https://brew.sh/ to install Homebrew.\n");
        print("Run this script again after installing Homebrew.\n");
        throw("Installation terminated because Homebrew is not installed yet");
    end
else
    brew = "";
end

print("Exasim find and install the required packages...\n");

gcc, appinstall = findinstallexec("g++", "gcc", brew, 0);
if appinstall==1
    gcc = findinstallexec("g++", "gcc", brew, 0);
end

mpi, appinstall = findinstallexec("mpicxx", "openmpi", brew, 0);
if appinstall==1
    mpi, appinstall = findinstallexec("mpicxx", "openmpi", brew, 0);
end

nvcc, appinstall = findinstallexec("nvcc", "nvcc", brew, 10);
if appinstall==1
    nvcc, appinstall = findinstallexec("nvcc", "nvcc", brew, 10);
end

metis, appinstall = findinstallexec("mpmetis", "metis", brew, 0);
if appinstall==1
    metis, appinstall = findinstallexec("mpmetis", "metis", brew, 0);
end

#python = findinstallexec("python", "python", brew, 0);
gmsh, appinstall = findinstallexec("gmsh", "gmsh", brew, 0);
if appinstall==1
    gmsh, appinstall = findinstallexec("gmsh", "gmsh", brew, 0);
end

if Sys.isapple()
    #julia = findinstallexec("julia", "julia", brew, 1);
    paraview,appinstall = findinstallexec("paraview", "paraview", brew, 1);
    if appinstall==1
        paraview, appinstall = findinstallexec("paraview", "paraview", brew, 1);
    end
else
    #julia = findinstallexec("julia", "julia", brew, 0);
    paraview,appinstall = findinstallexec("paraview", "paraview", brew, 0);
    if appinstall==1
        paraview, appinstall = findinstallexec("paraview", "paraview", brew, 0);
    end
end

if Sys.isapple()
    print("Installing " * "openblas" * " via brew.");
    run(string2cmd(brew * " install openblas"));
    print("Installing " * "lapack" * " via brew.");
    run(string2cmd(brew * " install lapack"));
else
    print("Installing " * "blas" * " via apt.");
    run(string2cmd("sudo apt install libblas-dev"));
    print("Installing " * "lapack" * " via apt.");
    run(string2cmd("sudo apt install liblapack-dev"));
end

cdir = pwd(); ii = findlast("Exasim", cdir);
versiondir = cdir[1:ii[end]] * "/" * version * "/Julia";

filem = versiondir * "/Preprocessing/initializepde.jl";
text = read(filem, String);

q = """ " """;
q = q[2];

oldgcc = "pde.cpucompiler = " * q * "g++" * q;
p = q * gcc * q;
newgcc = "pde.cpucompiler = " * p;
newtext = replace(text, oldgcc => newgcc);

oldmpi = "pde.mpicompiler = " * q * "mpicxx" * q;
p = q * mpi * q;
newmpi = "pde.mpicompiler = " * p;
newtext = replace(newtext, oldmpi => newmpi);

oldnvcc = "pde.gpucompiler = " * q * "nvcc" * q;
p = q * nvcc * q;
newnvcc = "pde.gpucompiler = " * p;
newtext = replace(newtext, oldnvcc => newnvcc);

oldmpirun = "pde.mpirun = " * q * "mpirun" * q;
p = q * mpi * q;
p = replace(p, "mpicxx" => "mpirun");
newmpirun = "pde.mpirun = " * p;
newtext = replace(newtext, oldmpirun => newmpirun);

oldmetis = "pde.metis = " * q * "mpmetis" * q;
p = q * metis * q;
newmetis = "pde.metis = " * p;
newtext = replace(newtext, oldmetis => newmetis);

oldgmsh = "pde.gmsh = " * q * "gmsh" * q;
p = q * gmsh * q;
newgmsh = "pde.gmsh = " * p;
newtext = replace(newtext, oldgmsh => newgmsh);

oldparaview = "pde.paraview = " * q * "paraview" * q;
p = q * paraview * q;
newparaview = "pde.paraview = " * p;
newtext = replace(newtext, oldparaview => newparaview);

fid = open(filem, "w");
write(fid, newtext);
close(fid);
