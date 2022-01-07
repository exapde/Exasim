# import external modules
import numpy, os
from numpy import pi

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
pde['torder'] = 3;          # time-stepping order of accuracy
pde['nstage'] = 3;          # time-stepping number of stages
pde['dt'] = 0.05*numpy.ones(200);   # time step sizes
pde['soltime'] = numpy.arange(1,pde['dt'].size+1); # steps at which solution are collected
pde['visdt'] = pde['dt'][0]; # visualization timestep size

gam = 1.4;                      # specific heat ratio
Re = 1600;                      # Reynolds number
Pr = 0.71;                      # Prandtl number
Minf = 0.2;                     # Mach number
pde['physicsparam'] = [gam, Re, Pr, Minf];
pde['tau'] = numpy.array([2.0]); # DG stabilization parameter
pde['GMRESrestart']=40;
pde['linearsolvertol']=0.0001;
pde['linearsolveriter']=41;
pde['precMatrixType']=2;
pde['NLiter']=2;

# create a mesh of 16 by 16 by 16 hexes on a cube domain
mesh['p'], mesh['t'] = Mesh.cubemesh(16,16,16,1)[0:2];
mesh['p'] = 2*pi*mesh['p'];
# expressions for domain boundaries
mesh['boundaryexpr'] = [lambda p: (p[1,:] < 1e-3), lambda p: (p[0,:] > 2*pi-1e-3), lambda p: (p[1,:] > 2*pi-1e-3), lambda p: (p[0,:] < 1e-3), lambda p: (p[2,:] < 1e-3), lambda p: (p[2,:] > 2*pi-1e-3)];
mesh['boundarycondition'] = numpy.array([1, 1, 1, 1, 1, 1]); # Set boundary condition for each boundary
mesh['periodicexpr'] = [[2, lambda p: p[[1,2],:], 4, lambda p: p[[1,2],:]], [1, lambda p: p[[0,2],:], 3, lambda p: p[[0,2],:]], [5, lambda p: p[[0,1],:], 6, lambda p: p[[0,1],:]]];

# call exasim to generate and run C++ code to solve the PDE model
sol, pde, mesh  = Postprocessing.exasim(pde,mesh)[0:3];

# visualize the numerical solution of the PDE model using Paraview
pde['visscalars'] = ["density", 0, "energy", 4]; # list of scalar fields for visualization
pde['visvectors'] = ["momentum", [1, 2, 3]]; # list of vector fields for visualization
Postprocessing.vis(sol,pde,mesh); # visualize the numerical solution
print("Done!");
