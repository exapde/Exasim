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
pde['model'] = "ModelC";       # ModelC, ModelD, ModelW
pde['modelfile'] = "pdemodel"; # name of a file defining the PDE model

# Choose computing platform and set number of processors
#pde['platform'] = "gpu";   # choose this option if NVIDIA GPUs are available
pde['mpiprocs'] = 1;        # number of MPI processors

# Set discretization parameters, physical parameters, and solver parameters
pde['porder'] = 4;         # polynomial degree
pde['torder'] = 3;          # time-stepping order of accuracy
pde['nstage'] = 3;          # time-stepping number of stages
pde['dt'] = 0.05*numpy.ones(200);   # time step sizes
pde['soltime'] = numpy.arange(10,pde['dt'].size+1,10); # steps at which solution are collected
pde['visdt'] = pde['dt'][0]; # visualization timestep size

gam = 5/3;                      # specific heat ratio
pde['physicsparam'] = [gam];
pde['tau'] = numpy.array([4.0]); # DG stabilization parameter

# create a mesh of 16 by 16 quads on a square domain
mesh['p'], mesh['t'] = Mesh.SquareMesh(16,16,1)[0:2];
mesh['p'] = 10*mesh['p'] - 5;
# expressions for domain boundaries
mesh['boundaryexpr'] = [lambda p: (p[1,:] < -5+1e-3), lambda p: (p[0,:] > 5-1e-3), lambda p: (p[1,:] > 5-1e-3), lambda p: (p[0,:] < -5+1e-3)];
mesh['boundarycondition'] = numpy.array([1, 1, 1, 1]); # Set boundary condition for each boundary
mesh['periodicexpr'] = [[2, lambda p: p[1,:], 4, lambda p: p[1,:]], [1, lambda p: p[0,:], 3, lambda p: p[0,:]]];

# call exasim to generate and run C++ code to solve the PDE model
sol, pde, mesh  = Postprocessing.exasim(pde,mesh)[0:3];

