# import external modules
import numpy, os

# Add Exasim to Python search path
cdir = os.getcwd(); ii = cdir.find("Exasim");
exec(open(cdir[0:(ii+6)] + "/install/setpath.py").read());

# import internal modules
import Preprocessing, Postprocessing, Gencode, Mesh

# Create pde object and mesh object
pde,mesh = Preprocessing.initializeexasim();

# Define a PDE model: governing equations and boundary conditions
pde['model'] = "ModelC";       # ModelC, ModelD, ModelW
pde['modelfile'] = "pdemodel"; # name of a file defining the PDE model

# Choose computing platform and set number of processors
#pde['platform'] = "gpu";   # choose this option if NVIDIA GPUs are available
pde['mpiprocs'] = 4;        # number of MPI processors

# Set discretization parameters, physical parameters, and solver parameters
pde['porder'] = 4;         # polynomial degree
pde['torder'] = 3;          # time-stepping order of accuracy
pde['nstage'] = 3;          # time-stepping number of stages
pde['physicsparam'] = [1.0, 1.0];    # convective velocity
pde['tau'] = numpy.array([1.0]);            # DG stabilization parameter
pde['dt'] = 0.025*numpy.ones(400);   # time step sizes
pde['soltime'] = numpy.arange(1,pde['dt'].size+1); # steps at which solution are collected
pde['visdt'] = 0.025; # visualization timestep size

# create a mesh of 10 by 10 quads on a square domain
mesh['p'], mesh['t'] = Mesh.SquareMesh(20,20,1)[0:2];
mesh['p'] = mesh['p'] - 0.5;
# expressions for domain boundaries
mesh['boundaryexpr'] = [lambda p: (p[1,:] < -0.5+1e-3), lambda p: (p[0,:] > 0.5-1e-3), lambda p: (p[1,:] > 0.5-1e-3), lambda p: (p[0,:] < -0.5+1e-3)];
mesh['boundarycondition'] = numpy.array([1, 1, 1, 1]); # Set boundary condition for each boundary

# call exasim to generate and run C++ code to solve the PDE model
sol, pde, mesh  = Postprocessing.exasim(pde,mesh)[0:3];

# visualize the numerical solution of the PDE model using Paraview
pde['visscalars'] = ["temperature", 0]; # list of scalar fields for visualization
pde['visvectors'] = []; # list of vector fields for visualization
Postprocessing.vis(sol,pde,mesh); # visualize the numerical solution
print("Done!");


