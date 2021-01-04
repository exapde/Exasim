# create pde structure and mesh structure
pde[2], mesh[2] = Preprocessing.initializeexasim(version);

# Define PDE model: governing equations, initial solutions, and boundary conditions
pde[2].model = "ModelD";            # ModelC, ModelD, ModelW
pde[2].modelfile = "pdemodel2";     # name of a file defining the PDE model
include(pde[2].modelfile * ".jl");  # include the PDE model file
pde[2].modelnumber = 2;

# Set discretization parameters, physical parameters, and solver parameters
pde[2].porder = 3;                  # polynomial degree
pde[2].physicsparam = [1.0 0.0];    # thermal conductivity and boundary value
pde[2].tau = [1.0];                 # DG stabilization parameter

# Choose computing platform and set number of processors
#pde[2].platform = "gpu";           # choose this option if NVIDIA GPUs are available
pde[2].mpiprocs = 1;                # number of MPI processors

# set indices to obtain v from the solutions of the other PDE models 
# first column : model index
# second column: solution index
pde[2].vindx = [1 1];

# create a linear mesh for a square domain
mesh[2].p, mesh[2].t = Mesh.SquareMesh(8,8,1); # a mesh of 8 by 8 quadrilaterals
# expressions for disjoint boundaries
mesh[2].boundaryexpr = [p -> (p[2,:] .< 1e-3), p -> (p[1,:] .> 1-1e-3), p -> (p[2,:] .> 1-1e-3), p -> (p[1,:] .< 1e-3)];
mesh[2].boundarycondition = [1 1 1 1]; # Set boundary condition for each disjoint boundary

