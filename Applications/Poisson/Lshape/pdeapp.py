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
pde['physicsparam'] = numpy.array([1.0]);   # unit thermal conductivity
pde['tau'] = numpy.array([1.0]);            # DG stabilization parameter

# call Gmsh to generate a mesh on L-shaped domain, see lshape.geo for details
dim = 2; elemtype = 0;
mesh['p'], mesh['t'] = Mesh.gmshcall(pde, "lshape", dim, elemtype)[0:2];
# expressions for domain boundaries
mesh['boundaryexpr'] = [lambda p: (p[1,:] < 2)];
mesh['boundarycondition'] = numpy.array([1]); # Set boundary condition for each boundary

# call exasim to generate and run C++ code to solve the PDE model
sol, pde, mesh  = Postprocessing.exasim(pde,mesh)[0:3];

# visualize the numerical solution of the PDE model using Paraview
pde['visscalars'] = ["temperature", 0]; # list of scalar fields for visualization
pde['visvectors'] = ["temperature gradient", numpy.array([1, 2]).astype(int)]; # list of vector fields for visualization
Postprocessing.vis(sol,pde,mesh); # visualize the numerical solution
print("Done!");
