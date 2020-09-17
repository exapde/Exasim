# Specify an Exasim version to run
version = "Version0.1";

# import external modules
import numpy, os

# Add Exasim to Python search path
cdir = os.getcwd(); ii = cdir.find("Exasim");
exec(open(cdir[0:(ii+6)] + "/Installation/setpath.py").read());

# import internal modules
import Preprocessing, Postprocessing, Gencode, Mesh

# Create pde object and mesh object
pde,mesh = Preprocessing.initializeexasim(version);

# Define a PDE model: governing equations and boundary conditions
pde['model'] = "ModelW";       # ModelC, ModelD, ModelW
pde['modelfile'] = "pdemodel"; # name of a file defining the PDE model

# Choose computing platform and set number of processors
#pde['platform'] = "gpu";   # choose this option if NVIDIA GPUs are available
pde['mpiprocs'] = 1;        # number of MPI processors

# Set discretization parameters, physical parameters, and solver parameters
pde['porder'] = 3;         # polynomial degree
pde['physicsparam'] = numpy.array([1.0]);   # unit thermal conductivity
pde['tau'] = numpy.array([1.0]);            # DG stabilization parameter

# read a grid from a file
mesh['p'], mesh['t'] = Mesh.readmesh("grid.bin",0);
mesh['t'] = mesh['t'] - 1;
# expressions for domain boundaries
mesh['boundaryexpr'] = [lambda p: (p[1,:] < -12+1e-3), lambda p: (p[0,:] > 12-1e-3), lambda p: (p[1,:] > 12-1e-3), lambda p: (p[0,:] < -12+1e-3), lambda p: (p[0,:] < 20)];
mesh['boundarycondition'] = numpy.array([1, 1, 1, 1, 2]); # Set boundary condition for each boundary
# expressions for curved boundaries
mesh['curvedboundary'] = numpy.array([0, 0, 0, 0, 1]);
mesh['curvedboundaryexpr'] = [lambda p: 0, lambda p: 0, lambda p: 0, lambda p: 0, lambda p: numpy.sqrt(p[0,:]*p[0,:]+p[1,:]*p[1,:])-1.0];

pde['elemtype'] = 0;
mesh['f'] = Preprocessing.facenumbering(mesh['p'],mesh['t'],pde['elemtype'],mesh['boundaryexpr'],mesh['periodicexpr'])[0];
mesh['dgnodes'] = Preprocessing.createdgnodes(mesh['p'],mesh['t'],mesh['f'],mesh['curvedboundary'],mesh['curvedboundaryexpr'],pde['porder']);

# # call exasim to generate and run C++ code to solve the PDE model
# sol, pde, mesh  = Postprocessing.exasim(pde,mesh)[0:3];
#
# # visualize the numerical solution of the PDE model using Paraview
# pde['visscalars'] = ["temperature", 0]; # list of scalar fields for visualization
# pde['visvectors'] = ["temperature gradient", numpy.array([1, 2]).astype(int)]; # list of vector fields for visualization
# Postprocessing.vis(sol,pde,mesh); # visualize the numerical solution
# print("Done!");
