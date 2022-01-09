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
pde.model = "ModelD";            # ModelC, ModelD, ModelW
include("pdemodel.jl");          # include the PDE model file

# Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;                  # polynomial degree
pde.physicsparam = [1.0 0.0];    # thermal conductivity and boundary value
pde.tau = [1.0];                 # DG stabilization parameter

# Choose computing platform and set number of processors
#pde.platform = "gpu";           # choose this option if NVIDIA GPUs are available
pde.mpiprocs = 2;                # number of MPI processors

# create a linear mesh of 8 by 8 by 8 quadrilaterals for a unit cube
mesh.p,mesh.t = Mesh.cubemesh(8,8,8,1);
# expressions for domain boundaries
mesh.boundaryexpr = [p -> (p[2,:] .< 1e-3), p -> (p[1,:] .> 1-1e-3), p -> (p[2,:] .> 1-1e-3), p -> (p[1,:] .< 1e-3), p -> (p[3,:] .< 1e-3), p -> (p[3,:] .> 1-1e-3)];
mesh.boundarycondition = [1 1 1 1 1 1]; # Set boundary condition for each disjoint boundary

# call exasim to generate and run C++ code to solve the PDE model
sol, pde, mesh, ~,~,~,~  = Postprocessing.exasim(pde,mesh);

# visualize the numerical solution of the PDE model using Paraview
#pde.paraview = "/Applications/ParaView-5.8.1.app/Contents/MacOS/paraview"; # Set the path to Paraview executable
pde.visscalars = ["temperature", 1];  # list of scalar fields for visualization
pde.visvectors = ["temperature gradient", [2, 3, 4]]; # list of vector fields for visualization
mesh.dgnodes = Postprocessing.vis(sol,pde,mesh); # visualize the numerical solution
x = mesh.dgnodes[:,1,:]; y = mesh.dgnodes[:,2,:]; z = mesh.dgnodes[:,3,:];
uexact = sin.(pi*x).*sin.(pi*y).*sin.(pi*z); # exact solution
uh = sol[:,1,:];  # numerical solution
maxerr = maximum(abs.(uh[:]-uexact[:]));
print("Maximum absolute error: $maxerr\n");
print("Done!");
