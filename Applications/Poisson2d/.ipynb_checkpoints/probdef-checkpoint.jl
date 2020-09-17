# specify an Exasim version to run
version = "Version0.1";

using Revise, DelimitedFiles, SymPy

# Add Exasim to Julia search path
cdir = pwd();
ii = findlast("Exasim", cdir);
versiondir = cdir[1:ii[end]] * "/" * version * "/Julia";
include(versiondir * "/setpath.jl");

using Structs, Master, Preprocessing, Mesh, Gencode

# create app struct
app = Structs.AppStruct();
app = Structs.InitializeAppStruct(app,version);

app.appname = "poi";   # application name
app.platform = "cpu";  # choose this option if NVIDIA GPUs are not available
#app.platform = "gpu";  # choose this option if NVIDIA GPUs are available
app.mpiprocs = 2;       # number of MPI ranks

# # #  Set CPU compiler # # #
# # #  use default g++ compiler on your system
app.cpucompiler = "g++";
# # #  or modify the below option to choose a specific CPU compiler
#  app.cpucompiler = "/opt/local/bin/g++-mp-10";

# # #  set CPU compiler options # # #
# # #  use the below options if blas/lapack libary is NOT in the system search path
#  app.cpuflags = "-O2 -ldl -lm -Wl,-rpath, /path/to/blaslapack -L/path/to/blaslapack -lblas -llapack";
# # #  use the below options if MKL libary is available on your system
#  app.cpuflags = "-O2 -ldl -lm -Wl,-rpath, /path/to/MKL -L/path/to/MKL -lmkl_intel_lp64 -lmkl_sequential -lmkl_core";

# # #  Set MPI compiler if app.mpiprocs > 1 # # #
#  app.mpicompiler = "mpicxx";
#  app.mpirun = "mpirun";
# # #  or modify the below options to choose a specific MPI compiler
 app.mpicompiler = "/opt/local/bin/mpicxx-openmpi-mp";
 app.mpirun = "/opt/local/bin/mpirun-openmpi-mp";

# # #  Set GPU compiler if app.platform = "gpu" # # #
# # #  use default nvcc compiler on your system
#  app.gpucompiler = "nvcc";
# # #  or modify the below option to choose a specific GPU compiler
#  app.gpucompiler = "/opt/local/bin/nvcc";

# # #  set GPU compiler options # # #
# # #  use the below options if CUDA libaries are NOT in the system search path
#  app.gpuflags = "-Wl,-rpath, /path/to/CUDAlibs -L/path/to/CUDAlibs -lcudart -lcublas";

# check if the above compilers exist on your system
Gencode.checkcompilers(app);

# Define PDE model, governing equations, and boundary conditions
app.pdemodel = 2;       # (u,q) type
app.porder = 3;         # polynomial degree

# Poisson equation with homogenous Dirichlet condition on a unit square
# q + \nabla u = 0 in \Omega
# \nabla dot flux = source in \Omega
# u = 0 on \partial Omega
# flux = 2*param*q
# source = 2*pi*pi*sin(pi*x)*sin(pi*y);
# udg = (u, qx, qy) with udg(1) = u, udg(2) = qx, udg(3) = qy
function flux(xdg, udg, odg, wdg, uinf, param)
    return param[1]*udg[2:end];
end
function source(xdg, udg, odg, wdg, uinf, param)
    return (2*pi*pi)*sin(pi*xdg[1])*sin(pi*xdg[2]);
end
function ubou(xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param)
    return 0.0;
end
function fbou(xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param)
    return param[1]*(udg[2]*nlg[1]+udg[3]*nlg[2])+tau[1]*(udg[1]-0.0);
end

app.Flux = "flux";      # name of the function defining PDE fluxes
app.Source = "source";  # name of the function defining source term
app.Fbou = "fbou";      # name of the function defining boundary flux
app.Ubou = "ubou";      # name of the function defining boundary value for the solution
app.dt = reshape([0.0],1,1);             # steady-state problem
app.boundaryconditions = [1 1 1 1]; # Set boundary condition for each boundary
app.physicsparam = reshape([1.0],1,1);   # unit thermal conductivity
app.tau = reshape([1.0],1,1);            # stabilization parameter
app.ncu = 1;

m = 8; n = 8; # m by n grid of the square
elemtype = 1; # quad elements
mesh = Structs.MeshStruct();
mesh.p,mesh.t = Mesh.SquareMesh(m,n,elemtype);
# expressions for domain boundaries
mesh.bndexpr = [p -> (p[2,:] .< 1e-3), p -> (p[1,:] .> 1-1e-3), p -> (p[2,:] .> 1-1e-3), p -> (p[1,:] .< 1e-3)];
# experssions for periodic boundaries
mesh.periodicexpr = [];

# sol struct
sol = Structs.SolStruct();
sol.dgnodes = Preprocessing.createnodes(mesh.p,mesh.t,app.porder);
sol.UDG = zeros(FloatP, size(sol.dgnodes,1), 3, size(sol.dgnodes,3));

# run preprocessing to save input data into binary files
app, master, dmd = Preprocessing.preprocessing(app,mesh.p,mesh.t,sol.dgnodes,sol.UDG,[],[],mesh.bndexpr,mesh.periodicexpr);

# generate source codes and store them in app
Gencode.gencode(app);

# compile source codes to build an executable file and store it in app
compilerstr = Gencode.compilecode(app);

# run executable file to compute solution and store it in dataout
runstr = Gencode.runcode(app);
