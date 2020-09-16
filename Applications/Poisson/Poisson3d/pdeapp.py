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
pde['model'] = "ModelD";       # ModelC, ModelD, ModelW
pde['modelfile'] = "pdemodel"; # name of a file defining the PDE model

# Set discretization parameters, physical parameters, and solver parameters
pde['porder'] = 3;         # polynomial degree
pde['physicsparam'] = numpy.array([1.0, 0.0]);   # unit thermal conductivity
pde['tau'] = numpy.array([1.0]);            # DG stabilization parameter

# Choose computing platform and set number of processors
#pde['platform = "gpu";   # choose this option if NVIDIA GPUs are available
pde['mpiprocs'] = 2;        # number of MPI processors

# create a mesh of 8 by 8 by 8 hexes for a unit cube
mesh['p'], mesh['t'] = Mesh.cubemesh(8,8,8,1)[0:2];
# expressions for domain boundaries
mesh['boundaryexpr'] = [lambda p: (p[1,:] < 1e-3), lambda p: (p[0,:] > 1-1e-3), lambda p: (p[1,:] > 1-1e-3), lambda p: (p[0,:] < 1e-3), lambda p: (p[2,:] < 1e-3), lambda p: (p[2,:] > 1-1e-3)];
mesh['boundarycondition'] = numpy.array([1, 1, 1, 1, 1, 1]); # Set boundary condition for each boundary

# call exasim to generate and run C++ code to solve the PDE model
sol, pde, mesh  = Postprocessing.exasim(pde,mesh)[0:3];

# visualize the numerical solution of the PDE model using Paraview
#pde['paraview'] = "/Applications/ParaView-5.8.1.app/Contents/MacOS/paraview"; # Set the path to Paraview executable
pde['visscalars'] = ["temperature", 0]; # list of scalar fields for visualization
pde['visvectors'] = ["temperature gradient", numpy.array([1, 2, 3]).astype(int)]; # list of vector fields for visualization
dgnodes = Postprocessing.vis(sol,pde,mesh); # visualize the numerical solution
x = dgnodes[:,0,:]; y = dgnodes[:,1,:]; z = dgnodes[:,2,:];
uexact = numpy.sin(numpy.pi*x)*numpy.sin(numpy.pi*y)*numpy.sin(numpy.pi*z); # exact solution
uh = sol[:,0,:];  # numerical solution
print("Maximum absolute error: %g\n" % max(abs(uh.flatten()-uexact.flatten())));
print("Done!");
