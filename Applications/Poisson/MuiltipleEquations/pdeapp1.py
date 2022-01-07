# Create pde object and mesh object
pde[0],mesh[0] = Preprocessing.initializeexasim(version);

# Define a PDE model: governing equations and boundary conditions
pde[0]['model'] = "ModelD";       # ModelC, ModelD, ModelW
pde[0]['modelfile'] = "pdemodel1"; # name of a file defining the PDE model
pde[0]['modelnumber'] = 1;

# Choose computing platform and set number of processors
#pde[0]['platform'] = "gpu";   # choose this option if NVIDIA GPUs are available
pde[0]['mpiprocs'] = 1;        # number of MPI processors

# Set discretization parameters, physical parameters, and solver parameters
pde[0]['porder'] = 3;         # polynomial degree
pde[0]['physicsparam'] = numpy.array([1.0]);   # unit thermal conductivity
pde[0]['tau'] = numpy.array([1.0]);            # DG stabilization parameter

# create a mesh of 8 by 8 quads on a square domain
mesh[0]['p'], mesh[0]['t'] = Mesh.SquareMesh(8,8,1)[0:2];
# expressions for domain boundaries
mesh[0]['boundaryexpr'] = [lambda p: (p[1,:] < 1e-3), lambda p: (p[0,:] > 1-1e-3), lambda p: (p[1,:] > 1-1e-3), lambda p: (p[0,:] < 1e-3)];
mesh[0]['boundarycondition'] = numpy.array([1, 1, 1, 1]); # Set boundary condition for each boundary





