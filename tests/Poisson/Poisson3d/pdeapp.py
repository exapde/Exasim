# import external modules
import numpy, os

# Add Exasim to Python search path
cdir = os.getcwd(); ii = cdir.find("Exasim");
exec(open(cdir[0:(ii+6)] + "/Installation/setpath.py").read());

# import internal modules
import Preprocessing, Postprocessing, Gencode, Mesh

# Create pde object and mesh object
pde,mesh = Preprocessing.initializeexasim();

# Define a PDE model: governing equations and boundary conditions
pde['model'] = "ModelD";       # ModelC, ModelD, ModelW
pde['modelfile'] = "pdemodel"; # name of a file defining the PDE model

# Set discretization parameters, physical parameters, and solver parameters
pde['porder'] = 3;         # polynomial degree
pde['physicsparam'] = numpy.array([1.0, 0.0]);   # unit thermal conductivity
pde['tau'] = numpy.array([1.0]);            # DG stabilization parameter

# Choose computing platform and set number of processors
#pde['platform'] = "gpu";   # choose this option if NVIDIA GPUs are available
pde['mpiprocs'] = 2;        # number of MPI processors

# create a mesh of 8 by 8 by 8 hexes for a unit cube
mesh['p'], mesh['t'] = Mesh.cubemesh(8,8,8,1)[0:2];
# expressions for domain boundaries
mesh['boundaryexpr'] = [lambda p: (p[1,:] < 1e-3), lambda p: (p[0,:] > 1-1e-3), lambda p: (p[1,:] > 1-1e-3), lambda p: (p[0,:] < 1e-3), lambda p: (p[2,:] < 1e-3), lambda p: (p[2,:] > 1-1e-3)];
mesh['boundarycondition'] = numpy.array([1, 1, 1, 1, 1, 1]); # Set boundary condition for each boundary

# call exasim to generate and run C++ code to solve the PDE model
sol, pde, mesh  = Postprocessing.exasim(pde,mesh)[0:3];

