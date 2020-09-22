# specify an Exasim version to run
version = "Version0.1";

# External packages
using Revise, DelimitedFiles, SymPy

# Add Exasim to Julia search path
cdir = pwd(); ii = findlast("Exasim", cdir);
include(cdir[1:ii[end]] * "/Installation/setpath.jl");

# Exasim packages
using Preprocessing, Mesh, Gencode, Postprocessing

# create pde structure and mesh structure
pde, mesh = Preprocessing.initializeexasim(version);

# Define PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelC";            # ModelC, ModelD, ModelW
include("pdemodel.jl");          # include the PDE model file

# Set discretization parameters, physical parameters, and solver parameters
pde.porder = 4;          # polynomial degree
pde.torder = 3;          # time-stepping order of accuracy
pde.nstage = 3;          # time-stepping number of stages
pde.dt = 0.05*ones(200);   # time step sizes
pde.soltime = collect(1:length(pde.dt)); # steps at which solution are collected
pde.visdt = 0.05; # visualization timestep size

gam = 1.4;                      # specific heat ratio
M_ref = sqrt(1/gam);            # Mach number
pde.physicsparam = [gam M_ref];
pde.tau = [1+1/M_ref];            # DG stabilization parameter

# Choose computing platform and set number of processors
#pde.platform = "gpu";           # choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;                # number of MPI processors

# create a linear mesh for a square domain
mesh.p, mesh.t = Mesh.SquareMesh(10,10,1); # a mesh of 8 by 8 quadrilaterals
mesh.p = 10.0*mesh.p .- 5.0;
# expressions for disjoint boundaries
mesh.boundaryexpr = [p -> (p[2,:] .< -5+1e-3), p -> (p[1,:] .> 5-1e-3), p -> (p[2,:] .> 5-1e-3), p -> (p[1,:] .< -5+1e-3)];
mesh.boundarycondition = [1 1 1 1]; # Set boundary condition for each disjoint boundary
mesh.periodicexpr = [2 p->p[2,:] 4 p->p[2,:]; 1 p->p[1,:] 3 p->p[1,:]];

# call exasim to generate and run C++ code to solve the PDE model
sol, pde, mesh,~,~,~,~  = Postprocessing.exasim(pde,mesh);

# visualize the numerical solution of the PDE model using Paraview
pde.visscalars = ["density", 1, "energy", 4];  # list of scalar fields for visualization
pde.visvectors = ["momentum", [2, 3]]; # list of vector fields for visualization
Postprocessing.vis(sol,pde,mesh); # visualize the numerical solution
print("Done!");
