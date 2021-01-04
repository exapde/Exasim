# create pde structure and mesh structure
pde[1], mesh[1] = Preprocessing.initializeexasim(version);

# Define PDE model: governing equations, initial solutions, and boundary conditions
pde[1].model = "ModelD";            # ModelC, ModelD, ModelW
pde[1].modelfile = "pdemodel1";    # name of a file defining the PDE model
include(pde[1].modelfile * ".jl"); # include the PDE model file
pde[1].modelnumber = 1;

# Set discretization parameters, physical parameters, and solver parameters
pde[1].porder = 3;                  # polynomial degree
pde[1].physicsparam = [1.0 0.0];    # thermal conductivity and boundary value
pde[1].tau = [1.0];                 # DG stabilization parameter

# Choose computing platform and set number of processors
#pde[1].platform = "gpu";           # choose this option if NVIDIA GPUs are available
pde[1].mpiprocs = 1;                # number of MPI processors

# create a linear mesh for a square domain
mesh[1].p, mesh[1].t = Mesh.SquareMesh(8,8,1); # a mesh of 8 by 8 quadrilaterals
# expressions for disjoint boundaries
mesh[1].boundaryexpr = [p -> (p[2,:] .< 1e-3), p -> (p[1,:] .> 1-1e-3), p -> (p[2,:] .> 1-1e-3), p -> (p[1,:] .< 1e-3)];
mesh[1].boundarycondition = [1 1 1 1]; # Set boundary condition for each disjoint boundary

