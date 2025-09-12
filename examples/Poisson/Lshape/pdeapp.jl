# External packages
using Revise, DelimitedFiles, SymPy

# Add Exasim to Julia search path
cdir = pwd(); ii = findlast("Exasim", cdir);
include(cdir[1:ii[end]] * "/install/setpath.jl");

# Exasim packages
using Preprocessing, Mesh, Gencode, Postprocessing

# create pde structure and mesh structure
pde, mesh = Preprocessing.initializeexasim();

# Define PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";            # ModelC, ModelD, ModelW
include("pdemodel.jl");          # include the PDE model file

# Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;                  # polynomial degree
pde.physicsparam = [1.0 0.0];    # thermal conductivity and boundary value
pde.tau = [1.0];                 # DG stabilization parameter

# Choose computing platform and set number of processors
#pde.platform = "gpu";           # choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;                # number of MPI processors

# call Gmsh to generate a mesh on L-shaped domain, see lshape.geo for details
dim = 2; elemtype = 0;
mesh.p, mesh.t = Mesh.gmshcall(pde, "lshape", dim, elemtype);
# expressions for disjoint boundaries
mesh.boundaryexpr = [p -> (p[2,:] .< 2)];
mesh.boundarycondition = reshape([1],1,1); # Set boundary condition for each disjoint boundary

# call exasim to generate and run C++ code to solve the PDE model
sol, pde, mesh,~,~,~,~  = Postprocessing.exasim(pde,mesh);

# visualize the numerical solution of the PDE model using Paraview
pde.visscalars = ["temperature", 1];  # list of scalar fields for visualization
pde.visvectors = ["temperature gradient", [2, 3]]; # list of vector fields for visualization
Postprocessing.vis(sol,pde,mesh); # visualize the numerical solution
print("Done!");
