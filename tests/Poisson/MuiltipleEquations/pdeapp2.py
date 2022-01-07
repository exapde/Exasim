# Create pde object and mesh object
pde[1],mesh[1] = Preprocessing.initializeexasim();

# Define a PDE model: governing equations and boundary conditions
pde[1]['model'] = "ModelD";       # ModelC, ModelD, ModelW
pde[1]['modelfile'] = "pdemodel2"; # name of a file defining the PDE model
pde[1]['modelnumber'] = 2;

# Choose computing platform and set number of processors
#pde[1]['platform'] = "gpu";   # choose this option if NVIDIA GPUs are available
pde[1]['mpiprocs'] = 1;        # number of MPI processors

# Set discretization parameters, physical parameters, and solver parameters
pde[1]['porder'] = 3;         # polynomial degree
pde[1]['physicsparam'] = numpy.array([1.0]);   # unit thermal conductivity
pde[1]['tau'] = numpy.array([1.0]);            # DG stabilization parameter

# set indices to obtain v from the solutions of the other PDE models 
# first column : model index
# second column: solution index
pde[1]['vindx'] = numpy.array([[1, 1]]);

# create a mesh of 8 by 8 quads on a square domain
mesh[1]['p'], mesh[1]['t'] = Mesh.SquareMesh(8,8,1)[0:2];
# expressions for domain boundaries
mesh[1]['boundaryexpr'] = [lambda p: (p[1,:] < 1e-3), lambda p: (p[0,:] > 1-1e-3), lambda p: (p[1,:] > 1-1e-3), lambda p: (p[0,:] < 1e-3)];
mesh[1]['boundarycondition'] = numpy.array([1, 1, 1, 1]); # Set boundary condition for each boundary





