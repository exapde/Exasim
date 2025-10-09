version = "src";

if ismac
    [brewstatus0,~] = system("which brew");
    [brewstatus1,~] = system("which /usr/bin/brew");
    [brewstatus2,~] = system("which /usr/local/bin/brew");
    if brewstatus0==0
        brew = "brew";
    elseif brewstatus1==0
        brew = "/usr/bin/brew";
    elseif brewstatus2==0
        brew = "/usr/local/bin/brew";
    else
        disp("Homebrew is not installed on your computer.")
        disp("Please visit https://brew.sh/ to install Homebrew.");
        disp("Run this script again after installing Homebrew.");
        disp("Installation terminated because Homebrew is not installed yet.");
        return;
    end    
else
    brew = "";
end

disp("Exasim find and install the required packages...");

[gcc, appinstall] = findinstallexec("g++", "gcc", brew, 0);
if appinstall==1
    gcc = findinstallexec("g++", "gcc", brew, 0);
end
    
[mpi, appinstall] = findinstallexec("mpicxx", "openmpi", brew, 0);
if appinstall==1
    mpi = findinstallexec("mpicxx", "openmpi", brew, 0);
end

[nvcc, appinstall] = findinstallexec("nvcc", "nvcc", brew, 10);
if appinstall==1
    nvcc = findinstallexec("nvcc", "nvcc", brew, 10);
end

[metis, appinstall] = findinstallexec("mpmetis", "metis", brew, 0);
if appinstall==1
    metis = findinstallexec("mpmetis", "metis", brew, 0);
end

%python = findinstallexec("python", "python", brew, 0);
[gmsh, appinstall] = findinstallexec("gmsh", "gmsh", brew, 0);
if appinstall==1
    gmsh = findinstallexec("gmsh", "gmsh", brew, 0);
end

if ismac            
    %julia = findinstallexec("julia", "julia", brew, 1);
    [paraview,appinstall] = findinstallexec("paraview", "paraview", brew, 1);
    if appinstall==1
        paraview = findinstallexec("paraview", "paraview", brew, 1);
    end    
else    
    %julia = findinstallexec("julia", "julia", brew, 0);
    [paraview,appinstall] = findinstallexec("paraview", "paraview", brew, 0);
    if appinstall==1
        paraview = findinstallexec("paraview", "paraview", brew, 0);
    end    
end

if ismac
    disp("Installing " + "openblas" + " via brew.");
    system(brew + " install openblas");                
    disp("Installing " + "lapack" + " via brew.");
    system(brew + " install lapack");                    
else
    disp("Installing " + "blas" + " via apt.");
    system("sudo apt install libblas-dev");                
    disp("Installing " + "lapack" + " via apt.");
    system("sudo apt install liblapack-dev");                          
end    

cdir = pwd(); ii = strfind(cdir, "Exasim");
versiondir = cdir(1:(ii+5)) + "/"  + version + "/Matlab";
filem = versiondir + "/Preprocessing/initializepde.m";
text = fileread(filem);

oldgcc = 'pde.cpucompiler = "g++"';
p = '"' + gcc + '"';
newgcc = 'pde.cpucompiler = ' + p;
newtext = strrep(text, oldgcc, newgcc);

oldmpi = 'pde.mpicompiler = "mpicxx"';
p = '"' + mpi + '"';
newmpi = 'pde.mpicompiler = ' + p;
newtext = strrep(newtext, oldmpi, newmpi);

oldnvcc = 'pde.gpucompiler = "nvcc"';
p = '"' + nvcc + '"';
newnvcc = 'pde.gpucompiler = ' + p;
newtext = strrep(newtext, oldnvcc, newnvcc);

oldmpirun = 'pde.mpirun = "mpirun"';
p = '"' + mpi + '"';
p = strrep(p, "mpicxx", "mpirun");
newmpirun = 'pde.mpirun = ' + p;
newtext = strrep(newtext, oldmpirun, newmpirun);

oldmetis = 'pde.metis = "mpmetis"';
p = '"' + metis + '"';
newmetis = 'pde.metis = ' + p;
newtext = strrep(newtext, oldmetis, newmetis);

oldgmsh = 'pde.gmsh = "gmsh"';
p = '"' + gmsh + '"';
newgmsh = 'pde.gmsh = ' + p;
newtext = strrep(newtext, oldgmsh, newgmsh);

oldparaview = 'pde.paraview = "paraview"';
p = '"' + paraview + '"';
newparaview = 'pde.paraview = ' + p;
newtext = strrep(newtext, oldparaview, newparaview);

fid = fopen(filem, "w");
fprintf(fid, char(newtext));
fclose(fid);





