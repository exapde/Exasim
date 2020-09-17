# specify an Exasim version to run
version = "Version0.1";

# External packages
using Revise, DelimitedFiles, SymPy

# Add Exasim to Julia search path
cdir = pwd(); ii = findlast("Exasim", cdir);
versiondir = cdir[1:ii[end]] * "/" * version * "/Julia";
include(versiondir * "/setpath.jl");

# Internal packages
using Structs, Master, Preprocessing, Mesh, Gencode

# create app struct
app = Structs.AppStruct();
app = Structs.InitializeAppStruct(app,version);

# Define PDE model: governing equations, initial solutions, and boundary conditions
app.appname = "poisson";       # application name
app.pdemodel = "ModelD";       # ModelC, ModelD, ModelW
app.pdemodelfile = "pdemodel.jl"; # name of a file defining the PDE model
include(app.pdemodelfile);

# Set discretization parameters, physical parameters, and solver parameters
app.porder = 3;         # polynomial degree
app.dt = reshape([0.0],1,1);             # steady-state problem
app.boundaryconditions = [1 1 1 1]; # Set boundary condition for each boundary
app.physicsparam = reshape([1.0],1,1);   # unit thermal conductivity
app.tau = reshape([1.0],1,1);            # DG stabilization parameter
app.ncu = 1;

app.platform = "cpu";    # choose this option if NVIDIA GPUs are not available
#app.platform = "gpu";   # choose this option if NVIDIA GPUs are available
app.mpiprocs = 1;        # number of MPI processors
app.cpucompiler = "g++"; # Clang/GNU C++ compiler
if app.mpiprocs>1        # set MPI compiler if app.mpiprocs > 1
    app.mpicompiler = "/opt/local/bin/mpicxx-openmpi-mp"; # path to a MPI compiler
    app.mpirun = "/opt/local/bin/mpirun-openmpi-mp";
end
if app.platform == "gpu"
    app.gpucompiler = "nvcc"; # Nvidia CUDA compiler (AMD GPU will be supported in the next version)
end
Gencode.checkcompilers(app);     # check if the compilers are available

# create a linear mesh
m = 8; n = 8; # m by n grid of the square
elemtype = 1; # quad elements
mesh = Structs.MeshStruct();
mesh.p,mesh.t = Mesh.SquareMesh(m,n,elemtype);
# expressions for domain boundaries
mesh.boundaryexpr = [p -> (p[2,:] .< 1e-3), p -> (p[1,:] .> 1-1e-3), p -> (p[2,:] .> 1-1e-3), p -> (p[1,:] .< 1e-3)];
# expressions for curved boundaries
mesh.curvedboundary = [0 0 0 0];
mesh.curvedboundaryexpr = [];
# experssions for periodic boundaries
mesh.periodicexpr = [];

# generate input files and store them in datain folder
app, master, dmd = Preprocessing.preprocessing(app,mesh);

# generate source codes and store them in app folder
Gencode.gencode(app);

# compile source codes to build an executable file and store it in app folder
compilerstr = Gencode.compilecode(app);

# run executable file to compute solution and store it in dataout folder
runstr = Gencode.runcode(app);
