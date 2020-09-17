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
pde.model = "ModelW";            # ModelC, ModelD, ModelW
include("pdemodel.jl");          # include the PDE model file

# Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;                  # polynomial degree
pde.physicsparam = [1.0 0.0];    # thermal conductivity and boundary value
pde.tau = [1.0];                 # DG stabilization parameter

# Choose computing platform and set number of processors
#pde.platform = "gpu";           # choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;                # number of MPI processors

# read a grid from a file
mesh.p, mesh.t = Mesh.readmesh("grid.bin",0);
# expressions for disjoint boundaries
mesh.boundaryexpr = [p -> (p[2,:] .< -12.0+1e-3), p -> (p[1,:] .> 12.0-1e-3), p -> (p[2,:] .> 12.0-1e-3), p -> (p[1,:] .< -12.0+1e-3), p -> (p[1,:] .< 20)];
mesh.boundarycondition = [1 1 1 1 2]; # Set boundary condition for each disjoint boundary
# expressions for curved boundaries
mesh.curvedboundary = [0 0 0 0 1];
mesh.curvedboundaryexpr = [p -> 0, p -> 0, p -> 0, p -> 0, p -> sqrt.(p[1,:].*p[1,:]+p[2,:].*p[2,:]).-1.0];

pde.elemtype = 0;
mesh.f,~,~ = Preprocessing.facenumbering(mesh.p,mesh.t,pde.elemtype,mesh.boundaryexpr,mesh.periodicexpr);
mesh.dgnodes = Preprocessing.createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);

# # call exasim to generate and run C++ code to solve the PDE model
# sol, pde, mesh,~,~,~,~  = Postprocessing.exasim(pde,mesh);
#
# # visualize the numerical solution of the PDE model using Paraview
# pde.visscalars = ["temperature", 1];  # list of scalar fields for visualization
# pde.visvectors = ["temperature gradient", [2, 3]]; # list of vector fields for visualization
# Postprocessing.vis(sol,pde,mesh); # visualize the numerical solution
# print("Done!");
