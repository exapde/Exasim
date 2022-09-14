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

# Choose computing platform and set number of processors
#pde['platform'] = "gpu";   # choose this option if NVIDIA GPUs are available
pde['mpiprocs'] = 1;        # number of MPI processors

# Set discretization parameters, physical parameters, and solver parameters
pde['porder'] = 3;         # polynomial degree
pde['physicsparam'] = [1.0, 2*numpy.pi];   # unit thermal conductivity
pde['tau'] = numpy.array([2*numpy.pi]);            # DG stabilization parameter


nDiv = 17
print(Mesh)
mesh['p'], mesh['t'] = Mesh.linemesh(nDiv)



# expressions for domain boundaries
mesh['boundaryexpr'] = [lambda p: (p[0,:] < 1e-3), lambda p: (p[0,:] > 1-1e-3)]
mesh['boundarycondition'] = numpy.array([1, 1]); # Set boundary condition for each boundary

# call exasim to generate and run C++ code to solve the PDE model
sol, pde, mesh  = Postprocessing.exasim(pde,mesh)[0:3]

print("Done!")



