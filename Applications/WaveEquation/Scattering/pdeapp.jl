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
pde.model = "ModelW";            # ModelC, ModelD, ModelW
include("pdemodel.jl");          # include the PDE model file

# Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;                  # polynomial degree
pde.torder = 3;                  # time-stepping order of accuracy
pde.nstage = 2;                  # time-stepping number of stages
pde.physicsparam = [1.0 10.0 0.0 1.0]; # relative permitivity, wave number vector, wave speed
pde.tau = [1.0];                 # DG stabilization parameter
pde.dt = 0.1*ones(100);   # time step sizes
pde.soltime = collect(1:length(pde.dt)); # steps at which solution are collected
pde.visdt = 0.1; # visualization timestep size
pde.GMRESrestart=40;
pde.NLiter=3;

# Choose computing platform and set number of processors
#pde.platform = "gpu";           # choose this option if NVIDIA GPUs are available
pde.mpiprocs = 4;                # number of MPI processors

# read a grid from a file
mesh.p, mesh.t = Mesh.readmesh("grid.bin",0);
# expressions for disjoint boundaries
mesh.boundaryexpr = [p -> (p[2,:] .< -12.0+1e-3), p -> (p[1,:] .> 12.0-1e-3), p -> (p[2,:] .> 12.0-1e-3), p -> (p[1,:] .< -12.0+1e-3), p -> (p[1,:] .< 20)];
mesh.boundarycondition = [1 1 1 1 2]; # Set boundary condition for each disjoint boundary
# expressions for curved boundaries
mesh.curvedboundary = [0 0 0 0 1];
mesh.curvedboundaryexpr = [p -> 0, p -> 0, p -> 0, p -> 0, p -> sqrt.(p[1,:].*p[1,:]+p[2,:].*p[2,:]).-1.0];

# call exasim to generate and run C++ code to solve the PDE model
sol, pde, mesh,~,~,~,~  = Postprocessing.exasim(pde,mesh);

# visualize the numerical solution of the PDE model using Paraview
pde.visscalars = ["velocity", 1, "displacement", 4];  # list of scalar fields for visualization
pde.visvectors = ["displacement gradient", [2, 3]]; # list of vector fields for visualization
Postprocessing.vis(sol,pde,mesh); # visualize the numerical solution
print("Done!");
