# External packages
using Revise, DelimitedFiles, SymPy

# Add Exasim to Julia search path
cdir = pwd(); ii = findlast("Exasim", cdir);
include(cdir[1:ii[end]] * "/Installation/setpath.jl");

# Exasim packages
using Preprocessing, Mesh, Gencode, Postprocessing

# create pde structure and mesh structure
pde, mesh = Preprocessing.initializeexasim();

# Define PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelC";            # ModelC, ModelD, ModelW
include("pdemodel.jl");          # include the PDE model file

# Set discretization parameters, physical parameters, and solver parameters
pde.porder = 4;          # polynomial degree
pde.torder = 3;          # time-stepping order of accuracy
pde.nstage = 3;          # time-stepping number of stages
pde.physicsparam = [1 1];    # convective velocity
pde.tau = [1.0];               # DG stabilization parameter
pde.dt = 0.025*ones(200);   # time step sizes
pde.soltime = collect(10:10:length(pde.dt)); # steps at which solution are collected
pde.visdt = 0.025; # visualization timestep size

# Choose computing platform and set number of processors
#pde.platform = "gpu";           # choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;                # number of MPI processors

# create a linear mesh for a square domain
mesh.p, mesh.t = Mesh.SquareMesh(20,20,1); # a mesh of 8 by 8 quadrilaterals
mesh.p = mesh.p .- 0.5; 
# expressions for disjoint boundaries
mesh.boundaryexpr = [p -> (p[2,:] .< -0.5+1e-3), p -> (p[1,:] .> 0.5-1e-3), p -> (p[2,:] .> 0.5-1e-3), p -> (p[1,:] .< -0.5+1e-3)];
mesh.boundarycondition = [1 1 1 1]; # Set boundary condition for each disjoint boundary

# call exasim to generate and run C++ code to solve the PDE model
sol, pde, mesh,~,~,~,~  = Postprocessing.exasim(pde,mesh);

